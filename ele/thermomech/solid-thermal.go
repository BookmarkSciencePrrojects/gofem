// Copyright 2015 Dorival Pedroso and Jaro Hokkanen. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package thermomech

import (
	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/tsr"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gofem/mdl/thermomech"
	elesolid "github.com/cpmech/gofem/ele/solid"
	mdlsolid "github.com/cpmech/gofem/mdl/solid"
	"github.com/cpmech/gofem/shp"
	"github.com/cpmech/gosl/utl"
	"github.com/cpmech/gofem/ele"
)

// ElemUT implements the solid-thermal coupling
type SolidThermal struct {

	// basic data
	Cell      *inp.Cell      // the cell structure
	LbbCell   *inp.Cell      // if LBB==false, same as Cell; otherwise LbbCell is a new cell with less vertices
	Edat      *inp.ElemData  // element data; stored in allocator to be used in Connect
	X         [][]float64    // matrix of nodal coordinates [ndim][nnode]
	Ndim      int            // space dimension
	Nu        int            // total number of unknown displacements
	Nt        int            // total number of unknown temperatures
	Nc        int            // if mixed formulation is used, number of fl variables

	// variables for dynamics
	Cdam      float64        // coefficient for damping // TODO: read this value
	Gfcn      fun.Func       // gravity function

	// optional data
	UseB      bool           // use B matrix
	Thickness float64        // thickness
	Debug     bool           // debugging flag

	// integration points (integration points of displacement element)
	IpsElem   []shp.Ipoint   // integration points of element
	IpsFace   []shp.Ipoint   // integration points corresponding to faces

	// material models and internal variables (sld model)
	SldMdl       mdlsolid.Model      // material model
	SldMdlSmall  mdlsolid.Small      // model specialisation for small strains
	SldMdlLarge  mdlsolid.Large      // model specialisation for large deformations
	TrmMdl       *thermomech.Thermomech      // thermal material model

	// internal variables
	States    []*mdlsolid.State   // [nip] states
	StatesBkp []*mdlsolid.State   // [nip] backup states
	StatesAux []*mdlsolid.State   // [nip] auxiliary backup states

	// additional variables
	Umap   []int             // assembly map (location array/element equations)
	Tmap   []int             // assembly map (location array/element equations)
	Cmap    []int            // [nc] map of "fl" variables (mixed formulation)
	NatBcs []*ele.NaturalBc  // natural boundary conditions
	Emat   [][]float64       // [nvert][nip] extrapolator matrix; if AddToExt is called

	// local starred variables for displacements
	ζs    [][]float64        // [nip][ndim] t2 star vars: ζ* = α1.u + α2.v + α3.a
	χs    [][]float64        // [nip][ndim] t2 star vars: χ* = α4.u + α5.v + α6.a
	divχs []float64          // [nip] divergent of χs (for coupled sims)

	// scratchpad. computed @ each ip (displacement element)
	grav []float64           // [ndim] gravity vector
	us   []float64           // [ndim] displacements @ ip
	fi      []float64        // [nu] internal forces
	B       [][]float64      // [nsig][nu] B matrix for axisymetric case
	D       [][]float64      // [nsig][nsig] constitutive consistent tangent matrix

	// scratchpad. computed @ each ip (temperature element)
	tstar   []float64        // [nip] ustar* = β1 u + β2 dudt
	xip     []float64        // real coordinates of ip
	tval    float64          // u(t,x) scalar field @ ip
	gradu   []float64        // [ndim] ∇u(t,x): gradient of u @ ip
	wvec    []float64        // [ndim] w(t,x) vector @ ip
	tmp     []float64        // auxiliary vector
	Sfun    fun.Func         // s(x) function for temperature element

	// strains
	ε       []float64        // total (updated) strains
	Δε      []float64        // incremental strains leading to updated strains

	//Mixed form
	MixedForm  bool          // indicates if this element has mixed form BC's
	Vid2convId []int         // [nverts] maps local vertex id to index in Cmap
	ConvId2vid []int         // [nf] maps mixed form BC face variable id to local vertex id

	// Stiffness matrices (u-t coupling)
	Kuu        [][]float64
	Kut        [][]float64
	Ktu        [][]float64
	Ktt        [][]float64

	//Additional stiffness matrices for mixed form solution
	Ktc        [][]float64   // [nt][nc]
	Kct        [][]float64   // [nc][nt]
	Kcc        [][]float64   // [nc][nc]
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	ele.SetInfoFunc("solid-thermal", func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *ele.Info {

		// number of nodes in element
		nverts := cell.Shp.Nverts
		if nverts < 0 {
			return nil // fail
		}

		// set DOFS and other information
		var info ele.Info
		ykeys := []string{"ux", "uy", "temp"}
		if sim.Ndim == 3 {
			ykeys = []string{"ux", "uy", "uz", "temp"}
		}
		info.Dofs = make([][]string, nverts)
		for m := 0; m < nverts; m++ {
			info.Dofs[m] = ykeys
			if m > cell.GetNverts(!sim.Data.NoLBB) - 1{
				info.Dofs[m] = info.Dofs[m][:len(info.Dofs[m])-1]
			}
		}

		info.Y2F = map[string]string{"ux": "fx", "uy": "fy", "uz": "fz", "temp" : "q"}

		// vertices on convective faces
		if len(cell.FaceBcs) > 0 {
			lverts := cell.FaceBcs.GetVerts("qcm")
			for _, m := range lverts {
				if m < cell.GetNverts(!sim.Data.NoLBB) { // avoid adding vertices of superelement (e.g. qua8 vertices in this qua4 cell)
					info.Dofs[m] = append(info.Dofs[m], "fl")
				}
			}
			if len(lverts) > 0 {
				ykeys = append(ykeys, "fl")
				info.Y2F["fl"] = "nil"
			}
		}

		info.T1vars = ykeys[sim.Ndim:]
		info.T2vars = ykeys[:sim.Ndim]

		// number of internal values to be extrapolated
		if cell.Extrap {
			info.Nextrap = 2 * sim.Ndim // nsig
		}

		return &info
	})

	// element allocator
	ele.SetAllocator("solid-thermal", func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) ele.Element {

		// basic data
		var o SolidThermal

		o.Cell = cell
		o.LbbCell = o.Cell

		if !sim.Data.NoLBB {
			o.LbbCell = o.Cell.GetSimilar(true)
		}

		o.Edat = edat
		o.Ndim = sim.Ndim
		o.Nu = o.Ndim * cell.Shp.Nverts
		o.Nt = o.LbbCell.Shp.Nverts
		o.X = x

		// vertices on convective faces
		var cverts []int
		if len(cell.FaceBcs) > 0 {
			lverts := cell.FaceBcs.GetVerts("qcm")
			for _, m := range lverts {
				if m < o.LbbCell.Shp.Nverts { // avoid adding vertices of superelement (e.g. qua8 vertices in this qua4 cell)
					cverts = append(cverts, m)
				}
			}
		}

		o.Nc = len(cverts)

		// parse flags
		o.UseB, o.Debug, o.Thickness = elesolid.GetSolidFlags(sim.Data.Axisym, sim.Data.Pstress, edat.Extra)

		// integration points
		var err error
		o.IpsElem, o.IpsFace, err = o.Cell.Shp.GetIps(edat.Nip, edat.Nipf)
		if err != nil {
			chk.Panic("cannot allocate integration points of ut element with nip=%d and nipf=%d:\n%v", edat.Nip, edat.Nipf, err)
		}
		nip := len(o.IpsElem)

		// model
		mat := sim.MatModels.Get(edat.Mat)
		if mat == nil {
			chk.Panic("cannot find material %q for solid-thermal element {tag=%d, id=%d}\n", edat.Mat, cell.Tag, cell.Id)
		}
		o.SldMdl = mat.Sld
		o.TrmMdl = mat.Trm.(*thermomech.Thermomech)

		// model specialisations
		switch m := o.SldMdl.(type) {
		case mdlsolid.Small:
			o.SldMdlSmall = m
		case mdlsolid.Large:
			o.SldMdlLarge = m
		default:
			chk.Panic("__internal_error__: 'u' element cannot determine the type of the material model")
		}

		// local starred variables
		o.ζs = la.MatAlloc(nip, o.Ndim)
		o.χs = la.MatAlloc(nip, o.Ndim)
		o.divχs = make([]float64, nip)

		// scratchpad. computed @ each ip
		nsig := 2 * o.Ndim
		o.grav = make([]float64, o.Ndim)
		o.us = make([]float64, o.Ndim)
		o.fi = make([]float64, o.Nu)
		o.D = la.MatAlloc(nsig, nsig)
		if o.UseB {
			o.B = la.MatAlloc(nsig, o.Nu)
		}

		// strains
		o.ε = make([]float64, nsig)
		o.Δε = make([]float64, nsig)

		// surface loads (natural boundary conditions)
		for _, fc := range cell.FaceBcs {
			o.NatBcs = append(o.NatBcs, &ele.NaturalBc{fc.Cond, fc.FaceId, fc.Func, fc.Extra})
		}

		// scratchpad
		o.tstar = make([]float64, nip)
		o.xip = make([]float64, o.Ndim)
		o.gradu = make([]float64, o.Ndim)
		o.wvec = make([]float64, o.Ndim)
		o.tmp = make([]float64, o.Ndim)

		o.Kuu = la.MatAlloc(o.Nu, o.Nu)
		o.Kut = la.MatAlloc(o.Nu, o.Nt)
		o.Ktu = la.MatAlloc(o.Nt, o.Nu)
		o.Ktt = la.MatAlloc(o.Nt, o.Nt)

		o.MixedForm = o.Nc > 0
		if o.MixedForm {

			// vertices on convection face; numbering
			o.ConvId2vid = cverts
			o.Vid2convId = utl.IntVals(o.Nt, -1)
			o.Cmap = make([]int, o.Nc)
			for μ, m := range o.ConvId2vid {
				o.Vid2convId[m] = μ
			}

			// allocate coupling matrices
			o.Ktc = la.MatAlloc(o.Nt, o.Nc)
			o.Kct = la.MatAlloc(o.Nc, o.Nt)
			o.Kcc = la.MatAlloc(o.Nc, o.Nc)
		}

		// return new element
		return &o
	})
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o *SolidThermal) Id() int { return o.Cell.Id }

// SetEqs set equations
func (o *SolidThermal) SetEqs(eqs [][]int, mixedform_eqs []int) (err error) {

	// Displacement DOFs
	o.Umap = make([]int, o.Nu)
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		for i := 0; i < o.Ndim; i++ {
			r := i + m*o.Ndim
			o.Umap[r] = eqs[m][i]
		}
	}
	// Temp DOFs
	o.Tmap = make([]int, o.Nt)
	for m := 0; m < o.LbbCell.Shp.Nverts; m++ {
		o.Tmap[m] = eqs[m][o.Ndim]
	}
	if o.MixedForm {
		for i, m := range o.ConvId2vid {
			o.Cmap[i] = eqs[m][o.Ndim+1]
		}
	}

	return
}

// SetEleConds set element conditions
func (o *SolidThermal) SetEleConds(key string, f fun.Func, extra string) (err error) {
	o.Sfun = nil
	switch key {

	case "g" : // gravity
		o.Gfcn = f
	case "s" : // source
		o.Sfun = f
	}

	return
}

// InterpStarVars interpolates star variables to integration points
func (o *SolidThermal) InterpStarVars(sol *ele.Solution) (err error) {

	// for each integration point
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}

		// clear local variables
		o.divχs[idx] = 0
		for i := 0; i < o.Ndim; i++ {
			o.ζs[idx][i], o.χs[idx][i] = 0, 0
		}

		// u-variables
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G
		for m := 0; m < o.Cell.Shp.Nverts; m++ {
			for i := 0; i < o.Ndim; i++ {
				r := o.Umap[i+m*o.Ndim]
				o.ζs[idx][i] += S[m] * sol.Zet[r]
				o.χs[idx][i] += S[m] * sol.Chi[r]
				o.divχs[idx] += G[m][i] * sol.Chi[r]
			}
		}

		// interpolation functions and gradients
		err = o.LbbCell.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}

		// clear local variables
		o.tstar[idx] =0

		// interpolate starred variables
		for m := 0; m < o.Nt; m++ {
			o.tstar[idx] += o.LbbCell.Shp.S[m] * sol.Psi[o.Tmap[m]]
		}
	}
	return
}

// adds -R to global residual vector fb
func (o *SolidThermal) AddToRhs(fb []float64, sol *ele.Solution) (err error) {

	if o.UseB {
		la.VecFill(o.fi, 0)
	}

	// for each integration point
	α1 := sol.DynCfs.GetAlp1()
	α4 := sol.DynCfs.GetAlp4()
	β1 := sol.DynCfs.GetBet1()
	prms := o.SldMdl.GetPrms()
	E , _ := prms.GetValues([]string{"E"})
	ρ := o.SldMdl.GetRho()
	u_nverts := o.Cell.Shp.Nverts
	t_nverts := o.LbbCell.Shp.Nverts
	var coef, dudt, kval, sval float64
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		err = o.ipvars(idx, sol)
		if err != nil {
			return
		}
		coef = o.Cell.Shp.J * ip[3]
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G
		St := o.LbbCell.Shp.S
		Gt := o.LbbCell.Shp.G

		// u residuals
		σe := o.States[idx].Sig
		if o.UseB {
			radius := 1.0
			if sol.Axisym {
				radius = o.Cell.Shp.AxisymGetRadius(o.X)
				coef *= radius
			}
			elesolid.IpBmatrix(o.B, o.Ndim, u_nverts, G, radius, S, sol.Axisym)
			la.MatTrVecMulAdd(o.fi, coef, o.B, σe) // fi += coef * tr(B) * σ
			for m := 0; m < u_nverts; m++ {
				for i := 0; i < o.Ndim; i++ {
					ru := o.Umap[i+m*o.Ndim]
					for n := 0; n < u_nverts; n++ {
						fb[ru] += coef * G[m][i] * S[n] * o.TrmMdl.Bcte[i] / E[0] * o.tval // coupling term Fu
					}
				}
			}
		} else {
			for m := 0; m < u_nverts; m++ {
				for i := 0; i < o.Ndim; i++ {
					ru := o.Umap[i+m*o.Ndim]
					fb[ru] -= coef * S[m] * (ρ*(α1*o.us[i]-o.ζs[idx][i]-o.grav[i]) + o.Cdam*(α4*o.us[i]-o.χs[idx][i])) // -RuBar
					for j := 0; j < o.Ndim; j++ {
						fb[ru] -= coef * tsr.M2T(σe, i, j) * G[m][j]
					}
					for n := 0; n < u_nverts; n++ {
						fb[ru] += coef * G[m][i] * S[n] * o.TrmMdl.Bcte[i] / E[0] * o.tval // coupling term Fu
					}
				}
			}
		}

		dudt = β1*o.tval - o.tstar[idx]
		kval = o.TrmMdl.Kval(o.tval)
		if o.Sfun != nil {
			sval = o.Sfun.F(sol.T, o.xip)
		}

		// compute wvec
		for i := 0; i < o.Ndim; i++ {
			o.wvec[i] = 0
			for j := 0; j < o.Ndim; j++ {
				o.wvec[i] -= kval * o.TrmMdl.Kcte[i][j] * o.gradu[j]
			}
		}

		// diffu residuals
		for m := 0; m < t_nverts; m++ {
			rt := o.Tmap[m]
			fb[rt] -= coef * St[m] * (ρ *dudt - sval) // - ftrs + fext/sval
			for i := 0; i < o.Ndim; i++ {
				fb[rt] += coef * Gt[m][i] * o.wvec[i] // + fint
				for n:=0; n < t_nverts; n++{
					fb[rt] += coef * Gt[m][i] * St[n] * o.TrmMdl.Bcte[i] / E[0] * (-α4 * o.us[i] + o.χs[idx][i]) // coupling term Ft
				}
			}
		}
	}

	// add fi term to fb, if using B matrix
	if o.UseB {
		for i, I := range o.Umap {
			fb[I] -= o.fi[i]
		}
	}

	// contribution from natural boundary conditions
	if len(o.NatBcs) > 0 {
		err = o.add_surfloads_to_rhs(fb, sol)
		if err != nil {
			return
		}
		err = o.add_natbcs_to_rhs(fb, sol)
		if err != nil {
			return
		}
	}

	return
}

// adds element K to global Jacobian matrix Kb
func (o *SolidThermal) AddToKb(Kb *la.Triplet, sol *ele.Solution, firstIt bool) (err error) {

	// zero K matrix
	la.MatFill(o.Kuu, 0)
	la.MatFill(o.Kut, 0)
	la.MatFill(o.Ktu, 0)
	la.MatFill(o.Ktt, 0)

	// for each integration point
	α1 := sol.DynCfs.GetAlp1()
	α4 := sol.DynCfs.GetAlp4()
	β1 := sol.DynCfs.GetBet1()
	prms := o.SldMdl.GetPrms()
	E , _ := prms.GetValues([]string{"E"})
	ρ := o.SldMdl.GetRho()
	u_nverts := o.Cell.Shp.Nverts
	t_nverts := o.LbbCell.Shp.Nverts
	var coef, dkdu, kval float64
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		err = o.ipvars(idx, sol)
		if err != nil {
			return
		}

		// check Jacobian
		if o.Cell.Shp.J < 0 {
			return chk.Err("ElemU: eid=%d: Jacobian is negative = %g\n", o.Id(), o.Cell.Shp.J)
		}

		// auxiliary
		coef = o.Cell.Shp.J * ip[3]
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G
		St := o.LbbCell.Shp.S
		Gt := o.LbbCell.Shp.G
		kval = o.TrmMdl.Kval(o.tval)
		dkdu = o.TrmMdl.DkDu(o.tval)

		// consistent tangent model matrix
		err = o.SldMdlSmall.CalcD(o.D, o.States[idx], firstIt)
		if err != nil {
			return
		}

		// add contribution to consistent tangent matrix
		if o.UseB {
			radius := 1.0
			if sol.Axisym {
				radius = o.Cell.Shp.AxisymGetRadius(o.X)
				coef *= radius
			}
			elesolid.IpBmatrix(o.B, o.Ndim, u_nverts, G, radius, S, sol.Axisym)
			la.MatTrMulAdd3(o.Kuu, coef, o.B, o.D, o.B) // K += coef * tr(B) * D * B
		} else {
			elesolid.IpAddToKt(o.Kuu, u_nverts, o.Ndim, coef, G, o.D)
		}

		// dynamic u term
		for m := 0; m < u_nverts; m++ {
			for i := 0; i < o.Ndim; i++ {
				r := i + m * o.Ndim
				for n := 0; n < u_nverts; n++ {
					c := i + n * o.Ndim
					o.Kuu[r][c] += coef * S[m] * S[n] * (ρ * α1 + o.Cdam * α4)
				}
			}
		}

		//diffu and coupling terms
		for m := 0; m < t_nverts; m++ {
			for i := 0; i < o.Ndim; i++ {
				o.tmp[i] = St[m]*dkdu*o.gradu[i] + kval*Gt[m][i]
			}
			for n := 0; n < u_nverts; n++ {
				for i := 0; i < o.Ndim; i++ {
					r := i + n * o.Ndim
					o.Kut[r][m] -= coef * Gt[m][i] * S[n] * o.TrmMdl.Bcte[i] / E[0]
					o.Ktu[m][r] += coef * G[n][i] * St[m] * o.TrmMdl.Bcte[i] / E[0] * α4
				}
			}
			for n := 0; n < t_nverts; n++ {
				o.Ktt[n][m] += coef * St[n] * St[m] * β1 * ρ
				for i := 0; i < o.Ndim; i++ {
					for j := 0; j < o.Ndim; j++ {
						o.Ktt[n][m] += coef * Gt[n][i] * o.TrmMdl.Kcte[i][j] * o.tmp[j]
					}
				}
			}
		}
	}

	// contribution from natural boundary conditions
	if len(o.NatBcs) > 0 {
		o.add_natbcs_to_jac(sol)
		if err != nil {
			return
		}
	}

	// add Ks to sparse matrix Kb

	for i, I := range o.Umap {
		for j, J := range o.Umap {
			Kb.Put(I, J, o.Kuu[i][j])
		}
	}
	for i, I := range o.Tmap {
		for j, J := range o.Umap {
			Kb.Put(I, J, o.Ktu[i][j])
			Kb.Put(J, I, o.Kut[j][i])
		}
	}
	for i, I := range o.Tmap {
		for j, J := range o.Tmap {
			Kb.Put(I, J, o.Ktt[i][j])
		}
	}

	if o.MixedForm {
		// Mixed form solution
		//    _             _
		//   |  Kuu Kut  0   |
		//   |  Ktu Ktt Ktc  |
		//   |_ 0   Kct Kcc _|
		//
		for i, I := range o.Tmap {
			for j, J := range o.Cmap {
				Kb.Put(I, J, o.Ktc[i][j])
				Kb.Put(J, I, o.Kct[j][i])
			}
		}
		for i, I := range o.Cmap {
			for j, J := range o.Cmap {
				Kb.Put(I, J, o.Kcc[i][j])
			}
		}
	}
	return
}

// Update perform (tangent) update
func (o *SolidThermal) Update(sol *ele.Solution) (err error) {
	// for each integration point
	nverts := o.Cell.Shp.Nverts
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G

		// compute strains
		if o.UseB {
			radius := 1.0
			if sol.Axisym {
				radius = o.Cell.Shp.AxisymGetRadius(o.X)
			}
			elesolid.IpBmatrix(o.B, o.Ndim, nverts, G, radius, S, sol.Axisym)
			elesolid.IpStrainsAndIncB(o.ε, o.Δε, 2*o.Ndim, o.Nu, o.B, sol.Y, sol.ΔY, o.Umap)
		} else {
			elesolid.IpStrainsAndInc(o.ε, o.Δε, nverts, o.Ndim, sol.Y, sol.ΔY, o.Umap, G)
		}

		// call model update => update stresses
		err = o.SldMdlSmall.Update(o.States[idx], o.ε, o.Δε, o.Id(), idx, sol.T)
		if err != nil {
			return chk.Err("Update failed (eid=%d, ip=%d)\nΔε=%v\n%v", o.Id(), idx, o.Δε, err)
		}
	}
	return
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// SetIniIvs sets initial ivs for given values in sol and ivs map
func (o *SolidThermal) SetIniIvs(sol *ele.Solution, ivs map[string][]float64) (err error) {
	// allocate slices of states
	nip := len(o.IpsElem)
	o.States = make([]*mdlsolid.State, nip)
	o.StatesBkp = make([]*mdlsolid.State, nip)
	o.StatesAux = make([]*mdlsolid.State, nip)

	// has specified stresses?
	_, has_sig := ivs["sx"]

	// for each integration point
	σ := make([]float64, 2*o.Ndim)
	for i := 0; i < nip; i++ {
		if has_sig {
			elesolid.Ivs2sigmas(σ, i, ivs)
		}
		o.States[i], err = o.SldMdl.InitIntVars(σ)
		if err != nil {
			return
		}
		o.StatesBkp[i] = o.States[i].GetCopy()
		o.StatesAux[i] = o.States[i].GetCopy()
	}
	return
}

// BackupIvs create copy of internal variables
func (o *SolidThermal) BackupIvs(aux bool) (err error) {
	if aux {
		for i, s := range o.StatesAux {
			s.Set(o.States[i])
		}
		return
	}
	for i, s := range o.StatesBkp {
		s.Set(o.States[i])
	}
	return
}

// RestoreIvs restore internal variables from copies
func (o *SolidThermal) RestoreIvs(aux bool) (err error) {
	if aux {
		for i, s := range o.States {
			s.Set(o.StatesAux[i])
		}
		return
	}
	for i, s := range o.States {
		s.Set(o.StatesBkp[i])
	}
	return
}

// Ureset fixes internal variables after u (displacements) have been zeroed
func (o *SolidThermal) Ureset(sol *ele.Solution) (err error) {
	for idx, _ := range o.IpsElem {
		if len(o.States[idx].F) > 0 {
			la.MatFill(o.States[idx].F, 0)
			la.MatFill(o.StatesBkp[idx].F, 0)
		}
	}
	return
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o *SolidThermal) Encode(enc utl.Encoder) (err error) {
	return enc.Encode(o.States)
}

// Decode decodes internal variables
func (o *SolidThermal) Decode(dec utl.Decoder) (err error) {
	err = dec.Decode(&o.States)
	if err != nil {
		return
	}
	return o.BackupIvs(false)
}

// OutIpCoords returns the coordinates of integration points
func (o *SolidThermal) OutIpCoords() (C [][]float64) {
	C = make([][]float64, len(o.IpsElem))
	for idx, ip := range o.IpsElem {
		C[idx] = o.Cell.Shp.IpRealCoords(o.X, ip)
	}
	return
}

// OutIpKeys returns the integration points' keys
func (o *SolidThermal) OutIpKeys() []string {
	keys := elesolid.StressKeys(o.Ndim)
	for i := 0; i < len(o.States[0].Alp); i++ {
		keys = append(keys, io.Sf("alp%d", i))
	}
	return keys
}

// OutIpVals returns the integration points' values corresponding to keys
func (o *SolidThermal) OutIpVals(M *ele.IpsMap, sol *ele.Solution) {
	nip := len(o.IpsElem)
	for i, key := range elesolid.StressKeys(o.Ndim) {
		for idx, _ := range o.IpsElem {
			M.Set(key, idx, nip, o.States[idx].Sig[i])
		}
	}
	for i := 0; i < len(o.States[0].Alp); i++ {
		key := io.Sf("alp%d", i)
		for idx, _ := range o.IpsElem {
			M.Set(key, idx, nip, o.States[idx].Alp[i])
		}
	}
}
// extra ////////////////////////////////////////////////////////////////////////////////////////////

// AddToExt extrapolates stresses at integration points to nodes
func (o *SolidThermal) AddToExt(sol *ele.Solution) (err error) {
	nverts := o.Cell.Shp.Nverts
	nsig := 2 * o.Ndim
	nip := len(o.IpsElem)
	if len(o.Emat) == 0 {
		o.Emat = la.MatAlloc(nverts, nip)
		err = o.Cell.Shp.Extrapolator(o.Emat, o.IpsElem)
		if err != nil {
			return
		}
	}
	for m := 0; m < nverts; m++ {
		vid := o.Cell.Verts[m]
		sol.Cnt[vid] += 1
		if len(sol.Ext[vid]) == 0 {
			sol.Ext[vid] = make([]float64, nsig)
		}
		for i := 0; i < nsig; i++ {
			for idx, _ := range o.IpsElem {
				σ := o.States[idx].Sig
				sol.Ext[vid][i] += o.Emat[m][idx] * σ[i]
			}
		}
	}
	return
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// ipvars computes current values @ integration points. idx == index of integration point
func (o *SolidThermal) ipvars(idx int, sol *ele.Solution) (err error) {

	// interpolation functions and gradients
	err = o.Cell.Shp.CalcAtIp(o.X, o.IpsElem[idx], true)
	if err != nil {
		return
	}
	// interpolation functions and gradients
	err = o.LbbCell.Shp.CalcAtIp(o.X, o.IpsElem[idx], true)
	if err != nil {
		return
	}

	// clear and recover t-variables @ ip
	o.tval = 0
	for i := 0; i < o.Ndim; i++ {
		o.us[i], o.gradu[i], o.xip[i] = 0, 0, 0
		for m := 0; m < o.Cell.Shp.Nverts; m++ {
			r := o.Umap[i+m*o.Ndim]
			o.us[i] += o.Cell.Shp.S[m] * sol.Y[r]
		}
	}

	// compute temp and its gradient @ ip by means of interpolating from nodes
	for m := 0; m < o.LbbCell.Shp.Nverts; m++ {
		r := o.Tmap[m]
		o.tval += o.LbbCell.Shp.S[m] * sol.Y[r]
		for i := 0; i < o.Ndim; i++ {
			o.gradu[i] += o.LbbCell.Shp.G[m][i] * sol.Y[r]
			o.xip[i] += o.LbbCell.Shp.S[m] * o.X[i][m]
		}
	}

	return
}

// surfloads_keys returns the keys that can be used to specify surface loads
func (o *SolidThermal) surfloads_keys() map[string]bool {
	return map[string]bool{"qn": true, "qn0": true, "aqn": true}
}

// add_surfloads_to_rhs adds surfaces loads to rhs
func (o *SolidThermal) add_surfloads_to_rhs(fb []float64, sol *ele.Solution) (err error) {

	// compute surface integral
	var res float64
	for _, nbc := range o.NatBcs {

		// function evaluation
		res = nbc.Fcn.F(sol.T, nil)

		// loop over ips of face
		for _, ipf := range o.IpsFace {

			// interpolation functions and gradients @ face
			iface := nbc.IdxFace
			err = o.Cell.Shp.CalcAtFaceIp(o.X, ipf, iface)
			if err != nil {
				return
			}
			Sf := o.Cell.Shp.Sf
			nvec := o.Cell.Shp.Fnvec

			// select natural boundary condition type
			switch nbc.Key {

			// distributed load
			case "qn", "qn0", "aqn":
				coef := ipf[3] * res * o.Thickness
				if sol.Axisym && nbc.Key == "aqn" {
					coef *= o.Cell.Shp.AxisymGetRadiusF(o.X, iface)
				}
				for j, m := range o.Cell.Shp.FaceLocalVerts[iface] {
					for i := 0; i < o.Ndim; i++ {
						r := o.Umap[i+m*o.Ndim]
						fb[r] += coef * Sf[j] * nvec[i] // +fe
					}
				}
			}
		}
	}
	return
}
// add_natbcs_to_rhs adds natural boundary conditions to rhs
func (o *SolidThermal) add_natbcs_to_rhs(fb []float64, sol *ele.Solution) (err error) {

	// compute surface integral
	var qb, temp0 float64
	for _, nbc := range o.NatBcs {

		// specified flux
		qb = nbc.Fcn.F(sol.T, nil)
		temp0 = nbc.Fcn.F(sol.T, nil)

		// loop over ips of face
		for _, ipf := range o.IpsFace {

			// interpolation functions and gradients @ face
			iface := nbc.IdxFace
			err = o.LbbCell.Shp.CalcAtFaceIp(o.X, ipf, iface)
			if err != nil {
				return
			}

			Sf := o.LbbCell.Shp.Sf
			Jf := la.VecNorm(o.LbbCell.Shp.Fnvec)
			coef := ipf[3] * Jf

			// select natural boundary condition type
			switch nbc.Key {

			// flux prescribed
			case "qb":
				for i, m := range o.LbbCell.Shp.FaceLocalVerts[iface] {
					fb[o.Tmap[m]] -= coef * qb * Sf[i]
				}
			case "qc":
				// clear temperature @ face ip
				tface := 0.0

				// compute tface @ face ip by means of interpolating from nodes
				for j, n := range o.LbbCell.Shp.FaceLocalVerts[iface] {
					tface += Sf[j] * sol.Y[o.Tmap[n]]
				}
				for i, m := range o.LbbCell.Shp.FaceLocalVerts[iface] {
					fb[o.Tmap[m]] -= coef * o.TrmMdl.H * (tface - temp0)  * Sf[i]
				}
			case "qcm":
				// clear temperature @ face ip
				tface := 0.0
				fl := 0.0

				// compute tface and fl @ face ip by means of interpolating from nodes
				for i, m := range o.LbbCell.Shp.FaceLocalVerts[iface] {
					μ := o.Vid2convId[m]
					tface += Sf[i] * sol.Y[o.Tmap[m]]
					fl += Sf[i] * sol.Y[o.Cmap[μ]]
				}
				for i, m := range o.LbbCell.Shp.FaceLocalVerts[iface] {
					μ := o.Vid2convId[m]
					fb[o.Tmap[m]] -= coef * o.TrmMdl.H * fl * Sf[i]
					fb[o.Cmap[μ]] -= coef * o.TrmMdl.H * (tface - fl - temp0)  * Sf[i]
				}
			}
		}
	}
	return
}
// add_natbcs_to_Ktt adds natural boundary conditions to Ktt
func (o *SolidThermal) add_natbcs_to_jac(sol *ele.Solution) (err error) {

	// compute surface integral
	for _, nbc := range o.NatBcs {

		// loop over ips of face
		for _, ipf := range o.IpsFace {

			// interpolation functions and gradients @ face
			iface := nbc.IdxFace
			err = o.LbbCell.Shp.CalcAtFaceIp(o.X, ipf, iface)
			if err != nil {
				return
			}

			Sf := o.LbbCell.Shp.Sf
			Jf := la.VecNorm(o.LbbCell.Shp.Fnvec)
			coef := ipf[3] * Jf

			// select natural boundary condition type
			switch nbc.Key {

			//Convection
			case "qc":
				for i, m := range o.LbbCell.Shp.FaceLocalVerts[iface] {
					for j, n := range o.LbbCell.Shp.FaceLocalVerts[iface] {
						o.Ktt[n][m] += coef * o.TrmMdl.H * Sf[i] * Sf[j]
					}
				}
			//Convection
			case "qcm":
				for i, m := range o.LbbCell.Shp.FaceLocalVerts[iface] {
					μ := o.Vid2convId[m]
					for j, n := range o.LbbCell.Shp.FaceLocalVerts[iface] {
						ν := o.Vid2convId[n]
						o.Ktc[m][ν] += coef * o.TrmMdl.H * Sf[i] * Sf[j]
						o.Kct[μ][n] += coef * o.TrmMdl.H * Sf[i] * Sf[j]
						o.Kcc[μ][ν] -= coef * o.TrmMdl.H * Sf[i] * Sf[j]
					}
				}
			}
		}
	}
	return
}
