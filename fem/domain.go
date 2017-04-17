// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
)

// Domain holds all Nodes and Elements active during a stage in addition to the Solution at nodes.
// Only elements in this processor are recorded here; however information from
// all cells might be recorded as well.
type Domain struct {

	// init: auxiliary variables
	Distr   bool            // distributed/parallel run
	Proc    int             // this processor number
	Verbose bool            // verbose
	ShowMsg bool            // show messages: if verbose==true and proc==0
	Sim     *inp.Simulation // [from FEM] input data
	Reg     *inp.Region     // region data
	Msh     *inp.Mesh       // mesh data
	LinSol  la.LinSol       // linear solver
	DynCfs  *ele.DynCoefs   // [from FEM] coefficients for dynamics/transient simulations

	// stage: nodes (active) and elements (active AND in this processor)
	Nodes  []*Node       // active nodes (for each stage). Note: indices in Nodes do NOT correpond to Ids => use Vid2node to access Nodes using Ids.
	Elems  []ele.Element // [procNcells] only active elements in this processor (for each stage)
	MyCids []int         // [procNcells] the ids of cells in this processor

	// stage: auxiliary maps for dofs and equation types
	F2Y      map[string]string // converts f-keys to y-keys; e.g.: "ux" => "fx"
	YandC    map[string]bool   // y and constraints keys; e.g. "ux", "pl", "H", "incsup", "rigid"
	Dof2Tnum map[string]int    // {t1,t2}-types: dof => t_number; e.g. "ux" => 2, "pl" => 1

	// stage: auxiliary maps for nodes and elements
	Vid2node   []*Node       // [nverts] VertexId => index in Nodes. Inactive vertices are 'nil'
	Cid2elem   []ele.Element // [ncells] CellId => index in Elems. Cells in other processors or inactive are 'nil'
	Cid2active []bool        // [ncells] CellId => whether cell is active or not in ANY processor

	// stage: subsets of elements
	ElemIntvars   []ele.WithIntVars    // elements with internal vars in this processor
	ElemIvsCon    []ele.WithIntVars    // elements with internal vars that are connectors
	ElemIvsNotCon []ele.WithIntVars    // elements with internal vars that are not connectors
	ElemConnect   []ele.Connector      // connector elements in this processor
	ElemExtrap    []ele.CanExtrapolate // elements with internal values to be extrapolated
	ElemFixedKM   []ele.WithFixedKM    // elements with fixed K,M matrices; to be recomputed if prms are changed

	// stage: coefficients and prescribed forces
	EssenBcs EssentialBcs // constraints (Lagrange multipliers)
	PtNatBcs PtNaturalBcs // point loads such as prescribed forces at nodes

	// stage: t1 and t2 variables
	T1eqs []int // first t-derivative variables; e.g.:  dp/dt vars (subset of ykeys)
	T2eqs []int // second t-derivative variables; e.g.: d²u/dt² vars (subset of ykeys)

	// stage: dimensions
	NnzKb int // number of nonzeros in Kb matrix
	Ny    int // total number of dofs, except λ
	Nlam  int // total number of Lagrange multipliers
	NnzA  int // number of nonzeros in A (constraints) matrix
	Nyb   int // total number of equations: ny + nλ

	// stage: solution and linear solver
	Sol      *ele.Solution // solution state
	Kb       *la.Triplet   // Jacobian == dRdy
	Fb       []float64     // residual == -fb
	Wb       []float64     // workspace
	InitLSol bool          // flag telling that linear solver needs to be initialised prior to any further call

	// for divergence control
	bkpSol *ele.Solution // backup solution
}

// Free frees memory
func (o *Domain) Free() {
	o.LinSol.Free()
	o.InitLSol = true // tell solver that lis has to be initialised before use
}

// NewDomains returns domains
func NewDomains(sim *inp.Simulation, dyncfs *ele.DynCoefs, proc, nproc int, distr, verb bool) (doms []*Domain) {
	doms = make([]*Domain, len(sim.Regions))
	for i, reg := range sim.Regions {
		doms[i] = new(Domain)
		doms[i].Distr = distr
		doms[i].Proc = proc
		doms[i].Verbose = verb
		doms[i].ShowMsg = verb && proc == 0
		doms[i].Sim = sim
		doms[i].Reg = reg
		doms[i].Msh = reg.Msh
		if distr {
			if nproc != len(reg.Msh.Part2cells) {
				chk.Panic("number of processors must be equal to the number of partitions defined in mesh file. %d != %d", nproc, len(reg.Msh.Part2cells))
			}
		}
		doms[i].LinSol = la.GetSolver(sim.LinSol.Name)
		doms[i].DynCfs = dyncfs
	}
	return
}

// SetStage set nodes, equation numbers and auxiliary data for given stage
func (o *Domain) SetStage(stgidx int) (err error) {

	// pointer to stage structure
	stg := o.Sim.Stages[stgidx]

	// backup state
	if stgidx > 0 {
		o.create_stage_copy()
		err = o.fix_inact_flags(stg.Activate, false)
		if err != nil {
			return
		}
		err = o.fix_inact_flags(stg.Deactivate, true)
		if err != nil {
			return
		}
	}

	// nodes (active) and elements (active AND in this processor)
	o.Nodes = make([]*Node, 0)
	o.Elems = make([]ele.Element, 0)
	o.MyCids = make([]int, 0)

	// auxiliary maps for dofs and equation types
	o.F2Y = make(map[string]string)
	o.YandC = GetIsEssenKeyMap()
	o.Dof2Tnum = make(map[string]int)

	// auxiliary maps for nodes and elements
	o.Vid2node = make([]*Node, len(o.Msh.Verts))
	o.Cid2elem = make([]ele.Element, len(o.Msh.Cells))
	o.Cid2active = make([]bool, len(o.Msh.Cells))

	// subsets of elements
	o.ElemConnect = make([]ele.Connector, 0)
	o.ElemIntvars = make([]ele.WithIntVars, 0)

	// allocate nodes and cells (active only) -------------------------------------------------------

	// for each cell
	var eq int // current equation number => total number of equations @ end of loop
	o.NnzKb = 0
	for _, cell := range o.Msh.Cells {

		// set cell's face boundary conditions
		err = cell.SetFaceConds(stg, o.Sim.Functions)
		if err != nil {
			return
		}

		// get element info
		info, inactive, err := ele.GetInfo(cell, o.Reg, o.Sim)
		if err != nil {
			return chk.Err("get element information failed:\n%v", err)
		}

		// skip inactive element
		if inactive || cell.Disabled {
			continue
		}
		o.Cid2active[cell.Id] = true

		// for non-joint elements, add new DOFs
		if !cell.IsJoint {
			chk.IntAssert(len(info.Dofs), len(cell.Verts))

			// store y and f information
			for ykey, fkey := range info.Y2F {
				o.F2Y[fkey] = ykey
				o.YandC[ykey] = true
				o.YandC[ykey+"_ini"] = true
			}

			// t1 and t2 equations
			for _, ykey := range info.T1vars {
				o.Dof2Tnum[ykey] = 1
			}
			for _, ykey := range info.T2vars {
				o.Dof2Tnum[ykey] = 2
			}

			// loop over nodes of this element
			var eNdof int // number of DOFs of this elmeent
			for j, v := range cell.Verts {

				// new or existent node
				var nod *Node
				if o.Vid2node[v] == nil {
					nod = NewNode(o.Msh.Verts[v])
					o.Vid2node[v] = nod
					o.Nodes = append(o.Nodes, nod)
				} else {
					nod = o.Vid2node[v]
				}

				// set DOFs and equation numbers
				for _, ukey := range info.Dofs[j] {
					eq = nod.AddDofAndEq(ukey, eq)
					eNdof += 1
				}
			}

			// number of non-zeros
			o.NnzKb += eNdof * eNdof
		}

		// allocate element
		mycell := cell.Part == o.Proc // cell belongs to this processor
		if mycell || !o.Distr {

			// new element
			ele, err := ele.New(cell, o.Reg, o.Sim)
			if err != nil {
				return chk.Err("new element failed:\n%v", err)
			}
			o.Cid2elem[cell.Id] = ele
			o.Elems = append(o.Elems, ele)
			o.MyCids = append(o.MyCids, ele.Id())

			// give equation numbers to new element
			eqs := make([][]int, len(cell.Verts))
			for j, v := range cell.Verts {
				for _, dof := range o.Vid2node[v].Dofs {
					eqs[j] = append(eqs[j], dof.Eq)
				}
			}
			err = ele.SetEqs(eqs, nil)
			if err != nil {
				return chk.Err("cannot set element equations:\n%v", err)
			}

			// subsets of elements
			o.add_element_to_subsets(info, ele)
		}
	}

	// connect elements (e.g. Joints)
	for _, e := range o.ElemConnect {
		nnz, err := e.Connect(o.Cid2elem, o.Msh.Cells[e.Id()])
		if err != nil {
			return chk.Err("cannot connect rod-joint elements with solids and rods:\n%v", err)
		}
		o.NnzKb += nnz
	}

	// element conditions, essential and natural boundary conditions --------------------------------

	// (re)set constraints and prescribed forces structures
	o.EssenBcs.Init(o.Sim.LiqMdl)
	o.PtNatBcs.Reset()

	// element conditions
	var fcn fun.TimeSpace
	for _, ec := range stg.EleConds {
		cells, ok := o.Msh.CellTag2cells[ec.Tag]
		if !ok {
			return chk.Err("cannot find cells with tag = %d to assign conditions", ec.Tag)
		}
		for _, cell := range cells {
			e := o.Cid2elem[cell.Id]
			if e != nil { // set conditions only for this processor's / active element
				for j, key := range ec.Keys {
					fcn, err = o.Sim.Functions.Get(ec.Funcs[j])
					if err != nil {
						return
					}
					e.SetEleConds(key, fcn, ec.Extra)
				}
			}
		}
	}

	// face essential boundary conditions
	for _, fc := range stg.FaceBcs {
		pairs, ok := o.Msh.FaceTag2cells[fc.Tag]
		if !ok {
			return chk.Err("cannot find faces with tag = %d to assign face boundary conditions", fc.Tag)
		}
		for _, pair := range pairs {
			cell := pair.C
			faceId := pair.Fid
			for _, bc := range cell.FaceBcs {
				if faceId == bc.FaceId {
					var nodes []*Node
					for _, vid := range bc.GlobalVerts {
						nodes = append(nodes, o.Vid2node[vid])
					}
					if o.YandC[bc.Cond] {
						err = o.EssenBcs.Set(bc.Cond, nodes, bc.Func, bc.Extra)
						if err != nil {
							return chk.Err("setting of essential (face) boundary conditions failed:\n%v", err)
						}
					}
				}
			}
		}
	}

	// vertex boundary conditions
	for _, nc := range stg.NodeBcs {
		verts, ok := o.Msh.VertTag2verts[nc.Tag]
		if !ok {
			return chk.Err("cannot find vertices with tag = %d to assign node boundary conditions", nc.Tag)
		}
		for _, v := range verts {
			if o.Vid2node[v.Id] != nil { // set BCs only for active nodes
				n := o.Vid2node[v.Id]
				for j, key := range nc.Keys {
					fcn, err = o.Sim.Functions.Get(nc.Funcs[j])
					if err != nil {
						return
					}
					if o.YandC[key] {
						o.EssenBcs.Set(key, []*Node{n}, fcn, nc.Extra)
					} else {
						o.PtNatBcs.Set(o.F2Y[key], n, fcn, nc.Extra)
					}
				}
			}
		}
	}

	// resize slices --------------------------------------------------------------------------------

	// t1 and t2 equations
	o.T1eqs = make([]int, 0)
	o.T2eqs = make([]int, 0)
	for _, nod := range o.Nodes {
		for _, dof := range nod.Dofs {
			switch o.Dof2Tnum[dof.Key] {
			case 1:
				o.T1eqs = append(o.T1eqs, dof.Eq)
			case 2:
				o.T2eqs = append(o.T2eqs, dof.Eq)
			default:
				chk.Panic("t1 and t2 equations are incorrectly set")
			}
		}
	}

	// size of arrays
	o.Ny = eq
	o.Nlam, o.NnzA = o.EssenBcs.Build(o.Ny)
	o.Nyb = o.Ny + o.Nlam

	// solution structure and linear solver
	o.Sol = new(ele.Solution)
	o.Sol.Steady = o.Sim.Data.Steady
	o.Sol.Axisym = o.Sim.Data.Axisym
	o.Sol.Pstress = o.Sim.Data.Pstress
	o.Sol.DynCfs = o.DynCfs

	// linear system and linear solver
	o.Kb = new(la.Triplet)
	o.Fb = make([]float64, o.Nyb)
	o.Wb = make([]float64, o.Nyb)
	o.Kb.Init(o.Nyb, o.Nyb, o.NnzKb+2*o.NnzA)
	o.InitLSol = true // tell solver that lis has to be initialised before use

	// allocate arrays
	o.Sol.Y = make([]float64, o.Ny)
	o.Sol.ΔY = make([]float64, o.Ny)
	o.Sol.L = make([]float64, o.Nlam)
	if !o.Sim.Data.Steady {
		o.Sol.Dydt = make([]float64, o.Ny)
		o.Sol.D2ydt2 = make([]float64, o.Ny)
		o.Sol.Psi = make([]float64, o.Ny)
		o.Sol.Zet = make([]float64, o.Ny)
		o.Sol.Chi = make([]float64, o.Ny)
	}

	// extrapolated values array
	o.Sol.Ext = make(map[int][]float64, 0)
	o.Sol.Cnt = make(map[int]int, 0)

	// message
	if o.ShowMsg {
		io.Pf(">> Steady=%v, Axisym=%v, Pstress=%v\n", o.Sol.Steady, o.Sol.Axisym, o.Sol.Pstress)
		io.Pf(">> Number of equations = %d\n", o.Ny)
		io.Pf(">> Number of Lagrange multipliers = %d\n", o.Nlam)
	}

	// success
	return
}

// SetIniVals sets/resets initial values (nodes and integration points)
func (o *Domain) SetIniVals(stgidx int, zeroSol bool) (err error) {

	// pointer to stage structure
	stg := o.Sim.Stages[stgidx]

	// clear solution vectors
	if zeroSol {
		o.Sol.Reset(o.Sim.Data.Steady)
	}

	// initialise internal variables
	if stg.IniPorous != nil {
		err = o.IniSetPorous(stg)
		if err != nil {
			return
		}
		if o.ShowMsg {
			io.Pf(">> Initial porous media state set\n")
		}
	} else if stg.IniStress != nil {
		err = o.IniSetStress(stg)
		if err != nil {
			return
		}
		if o.ShowMsg {
			io.Pf(">> Initial stress state set\n")
		}
	} else if stg.IniFcn != nil {
		err = o.IniSetFileFunc(stg)
		if err != nil {
			return
		}
		if o.ShowMsg {
			io.Pf(">> Initial state set by using function\n")
		}
	} else {
		for _, e := range o.ElemIntvars {
			e.SetIniIvs(o.Sol, nil)
		}
		if o.ShowMsg {
			io.Pf(">> Initial state set with default values\n")
		}
	}

	// import results from another set of files
	if stg.IniImport != nil {
		sum := new(Summary)
		err = sum.Read(stg.IniImport.Dir, stg.IniImport.Fnk, o.Sim.EncType)
		if err != nil {
			return chk.Err("cannot import state from %s/%s.sim:\n%v", stg.IniImport.Dir, stg.IniImport.Fnk, err)
		}
		err = o.Read(sum, len(sum.OutTimes)-1, o.Proc, false)
		if err != nil {
			return chk.Err("cannot load results into domain:\n%v", err)
		}
		if o.Ny != len(o.Sol.Y) {
			return chk.Err("length of primary variables vector imported is not equal to the one allocated. make sure the number of DOFs of the imported simulation matches this one. %d != %d", o.Ny, len(o.Sol.Y))
		}
		if stg.IniImport.ResetU {
			for _, ele := range o.ElemIntvars {
				err = ele.Ureset(o.Sol)
				if err != nil {
					return chk.Err("cannot run reset function of element after displacements are zeroed:\n%v", err)
				}
			}
			for _, nod := range o.Nodes {
				for _, ukey := range []string{"ux", "uy", "uz"} {
					eq := nod.GetEq(ukey)
					if eq >= 0 {
						o.Sol.Y[eq] = 0
						if len(o.Sol.Dydt) > 0 {
							o.Sol.Dydt[eq] = 0
							o.Sol.D2ydt2[eq] = 0
						}
					}
				}
			}
		}
		if o.ShowMsg {
			io.Pf(">> Initial state overridden by 'import' file\n")
		}
	}

	// set boundary conditions that depend on initial values
	o.EssenBcs.FixIniVals(o.Sol)

	// list boundary conditions
	if o.Sim.Data.ListBcs {
		io.Pf("%v", o.EssenBcs.List(stg.Control.Tf))
	}

	// make sure time is zero at the beginning of simulation
	o.Sol.T = 0
	return
}

// UpdateElems update elements after Solution has been updated
func (o *Domain) UpdateElems() (err error) {

	// update elements with internal values that aren't connectors
	for _, e := range o.ElemIvsNotCon {
		err = e.Update(o.Sol)
		if err != nil {
			break
		}
	}

	// compute extrapolated values @ nodes
	if len(o.ElemExtrap) > 0 {
		for vid, val := range o.Sol.Ext {
			la.VecFill(val, 0)
			o.Sol.Cnt[vid] = 0
		}
		for _, e := range o.ElemExtrap {
			err = e.AddToExt(o.Sol)
			if err != nil {
				break
			}
		}
		for vid, val := range o.Sol.Ext {
			count := float64(o.Sol.Cnt[vid])
			for i := 0; i < len(val); i++ {
				val[i] /= count
			}
		}
		// TODO: join Ext from multiple processors
	}

	// update elements with internal values that are connectors
	for _, e := range o.ElemIvsCon {
		err = e.Update(o.Sol)
		if err != nil {
			break
		}
	}
	return
}

// RecomputeKM recompute K and M matrices of elements with static matrices
func (o *Domain) RecomputeKM() {
	for _, e := range o.ElemFixedKM {
		e.Recompute(!o.Sim.Data.Steady)
	}
}

// auxiliary functions //////////////////////////////////////////////////////////////////////////////

// add_element_to_subsets adds an Elem to many subsets as it fits
func (o *Domain) add_element_to_subsets(info *ele.Info, element ele.Element) {

	// elements with internal variables
	if e, ok := element.(ele.WithIntVars); ok {
		o.ElemIntvars = append(o.ElemIntvars, e)
		if _, connector := element.(ele.Connector); connector {
			o.ElemIvsCon = append(o.ElemIvsCon, e)
		} else {
			o.ElemIvsNotCon = append(o.ElemIvsNotCon, e)
		}
	}

	// connector elements
	if e, ok := element.(ele.Connector); ok {
		o.ElemConnect = append(o.ElemConnect, e)
	}

	// subset of elements with data to be extrapolated
	if info.Nextrap > 0 {
		if e, ok := element.(ele.CanExtrapolate); ok {
			o.ElemExtrap = append(o.ElemExtrap, e)
		}
	}

	// subset of elements with fixed KM
	if e, ok := element.(ele.WithFixedKM); ok {
		o.ElemFixedKM = append(o.ElemFixedKM, e)
	}
}

// create_stage_copy creates a copy of current stage => to be used later when activating/deactivating elements
func (o *Domain) create_stage_copy() {
}

// set_act_deact_flags sets inactive flags for new active/inactive elements
func (o *Domain) fix_inact_flags(eids_or_tags []int, deactivate bool) (err error) {
	for _, tag := range eids_or_tags {
		if tag >= 0 { // this meahs that tag == cell.Id
			cell := o.Msh.Cells[tag]
			tag = cell.Tag
		}
		edat := o.Reg.Etag2data(tag)
		if edat == nil {
			return chk.Err("cannot get element's data with etag=%d", tag)
		}
		edat.Inact = deactivate
	}
	return
}

// backup saves a copy of solution
func (o *Domain) backup() {
	if o.bkpSol == nil {
		o.bkpSol = new(ele.Solution)
		o.bkpSol.Y = make([]float64, o.Ny)
		o.bkpSol.ΔY = make([]float64, o.Ny)
		o.bkpSol.L = make([]float64, o.Nlam)
		if !o.Sim.Data.Steady {
			o.bkpSol.Dydt = make([]float64, o.Ny)
			o.bkpSol.D2ydt2 = make([]float64, o.Ny)
			o.bkpSol.Psi = make([]float64, o.Ny)
			o.bkpSol.Zet = make([]float64, o.Ny)
			o.bkpSol.Chi = make([]float64, o.Ny)
		}
	}
	o.bkpSol.T = o.Sol.T
	copy(o.bkpSol.Y, o.Sol.Y)
	copy(o.bkpSol.ΔY, o.Sol.ΔY)
	copy(o.bkpSol.L, o.Sol.L)
	if !o.Sim.Data.Steady {
		copy(o.bkpSol.Dydt, o.Sol.Dydt)
		copy(o.bkpSol.D2ydt2, o.Sol.D2ydt2)
		copy(o.bkpSol.Psi, o.Sol.Psi)
		copy(o.bkpSol.Zet, o.Sol.Zet)
		copy(o.bkpSol.Chi, o.Sol.Chi)
	}
	for _, e := range o.ElemIntvars {
		e.BackupIvs(true)
	}
}

// restore restores solution
func (o *Domain) restore() {
	o.Sol.T = o.bkpSol.T
	copy(o.Sol.Y, o.bkpSol.Y)
	copy(o.Sol.ΔY, o.bkpSol.ΔY)
	copy(o.Sol.L, o.bkpSol.L)
	if !o.Sim.Data.Steady {
		copy(o.Sol.Dydt, o.bkpSol.Dydt)
		copy(o.Sol.D2ydt2, o.bkpSol.D2ydt2)
		copy(o.Sol.Psi, o.bkpSol.Psi)
		copy(o.Sol.Zet, o.bkpSol.Zet)
		copy(o.Sol.Chi, o.bkpSol.Chi)
	}
	for _, e := range o.ElemIntvars {
		e.RestoreIvs(true)
	}
}
