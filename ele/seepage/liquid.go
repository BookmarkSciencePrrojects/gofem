// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package seepage implements elements for solving seepage problems
package seepage

import (
	"math"

	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mdl/porous"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// Liquid implements an element for transient seepage analyses [1]
//  References:
//   [1] Pedroso DM (2015) A solution to transient seepage in unsaturated porous media.
//       Computer Methods in Applied Mechanics and Engineering, 285 791-816,
//       http://dx.doi.org/10.1016/j.cma.2014.12.009
type Liquid struct {

	// basic data
	Cell *inp.Cell   // the cell structure
	X    [][]float64 // matrix of nodal coordinates [ndim][nnode]
	Np   int         // total number of unknowns == number of vertices
	Ndim int         // space dimension

	// integration points
	IpsElem []shp.Ipoint // integration points of element
	IpsFace []shp.Ipoint // integration points corresponding to faces

	// material model
	Mdl *porous.Model // model

	// problem variables
	Pmap []int // assembly map (location array/element equations)

	// internal variables
	States    []*porous.State
	StatesBkp []*porous.State
	StatesAux []*porous.State

	// gravity
	Gfcn fun.Func // gravity function

	// natural boundary conditions
	NatBcs []*ele.NaturalBc // natural boundary conditions

	// flux boundary conditions (qb == \bar{q})
	Rhol_ex     []float64   // [nverts] ρl extrapolted to nodes => if has qb (flux)
	Drholdpl_ex [][]float64 // [nverts][nverts] Cpl extrapolted to nodes => if has qb (flux)
	Emat        [][]float64 // [nverts][nips] extrapolator matrix
	DoExtrap    bool        // do extrapolation of ρl and Cpl => for use with flux and seepage conditions

	// seepage face
	Nf         int         // number of fl variables
	HasSeep    bool        // indicates if this element has seepage faces
	Vid2seepId []int       // [nverts] maps local vertex id to index in Fmap
	SeepId2vid []int       // [nf] maps seepage face variable id to local vertex id
	Fmap       []int       // [nf] map of "fl" variables (seepage face)
	Macaulay   bool        // use discrete ramp function instead of smooth ramp
	BetRmp     float64     // coefficient for Sramp
	Kap        float64     // κ coefficient to normalise equation for seepage face modelling
	Hst        []bool      // [nf] set hydrostatic plmax
	Plmax      [][]float64 // [nf][nipsFace] specified plmax (not corrected by multiplier)

	// local starred variables
	PsiL []float64 // [nip] ψl* = β1.p + β2.dpdt

	// scratchpad. computed @ each ip
	Grav   []float64      // [ndim] gravity vector
	Pl     float64        // pl: liquid pressure
	GradPl []float64      // [ndim] ∇pl: gradient of liquid pressure
	Rhowl  []float64      // [ndim] ρl*wl: weighted liquid relative velocity
	Kpp    [][]float64    // [np][np] Kpp := dRpl/dpl consistent tangent matrix
	Kpf    [][]float64    // [np][nf] Kpf := dRpl/dfl consistent tangent matrix
	Kfp    [][]float64    // [nf][np] Kfp := dRfl/dpl consistent tangent matrix
	Kff    [][]float64    // [nf][nf] Kff := dRfl/dfl consistent tangent matrix
	LsVars *porous.LsVars // variable to hold results from CalcLs
	Tmp    []float64      // [ndim] temporary (auxiliary) vector
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	ele.SetInfoFunc("liquid", func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *ele.Info {

		// new info
		var info ele.Info

		// number of nodes in element
		nverts := cell.GetNverts(edat.Lbb)

		// solution variables
		ykeys := []string{"pl"}
		info.Dofs = make([][]string, nverts)
		for m := 0; m < nverts; m++ {
			info.Dofs[m] = ykeys
		}

		// maps
		info.Y2F = map[string]string{"pl": "ql"}

		// vertices on seepage faces
		if len(cell.FaceBcs) > 0 {
			lverts := cell.FaceBcs.GetVerts("seep")
			for _, m := range lverts {
				if m < nverts { // avoid adding vertices of superelement (e.g. qua8 vertices in this qua4 cell)
					info.Dofs[m] = append(info.Dofs[m], "fl")
				}
			}
			if len(lverts) > 0 {
				ykeys = append(ykeys, "fl")
				info.Y2F["fl"] = "nil"
			}
		}

		// t1 and t2 variables
		info.T1vars = ykeys
		return &info
	})

	// element allocator
	ele.SetAllocator("liquid", func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) ele.Element {

		// basic data
		var o Liquid
		o.Cell = cell
		o.X = x
		o.Np = o.Cell.Shp.Nverts
		o.Ndim = sim.Ndim

		// integration points
		var err error
		o.IpsElem, o.IpsFace, err = o.Cell.Shp.GetIps(edat.Nip, edat.Nipf)
		if err != nil {
			chk.Panic("cannot allocate integration points of p-element with nip=%d and nipf=%d:\n%v", edat.Nip, edat.Nipf, err)
		}
		nip := len(o.IpsElem)

		// model
		mat := sim.MatModels.Get(edat.Mat)
		if mat == nil {
			chk.Panic("cannot get model for p-element {tag=%d id=%d material=%q}:\n%v", cell.Tag, cell.Id, edat.Mat, err)
		}
		o.Mdl = mat.Por

		// local starred variables
		o.PsiL = make([]float64, nip)

		// scratchpad. computed @ each ip
		o.Grav = make([]float64, o.Ndim)
		o.GradPl = make([]float64, o.Ndim)
		o.Rhowl = make([]float64, o.Ndim)
		o.Tmp = make([]float64, o.Ndim)
		o.Kpp = la.MatAlloc(o.Np, o.Np)
		o.LsVars = new(porous.LsVars)

		// vertices on seepage faces
		var seepverts []int
		if len(cell.FaceBcs) > 0 {
			lverts := cell.FaceBcs.GetVerts("seep")
			for _, m := range lverts {
				if m < o.Np { // avoid adding vertices of superelement (e.g. qua8 vertices in this qua4 cell)
					seepverts = append(seepverts, m)
				}
			}
		}

		o.Nf = len(seepverts)
		o.HasSeep = o.Nf > 0
		if o.HasSeep {

			// vertices on seepage face; numbering
			o.SeepId2vid = seepverts
			o.Vid2seepId = utl.IntVals(o.Np, -1)
			o.Fmap = make([]int, o.Nf)
			for μ, m := range o.SeepId2vid {
				o.Vid2seepId[m] = μ
			}

			// flags
			o.Macaulay, o.BetRmp, o.Kap = GetSeepFaceFlags(edat.Extra)

			// allocate coupling matrices
			o.Kpf = la.MatAlloc(o.Np, o.Nf)
			o.Kfp = la.MatAlloc(o.Nf, o.Np)
			o.Kff = la.MatAlloc(o.Nf, o.Nf)
		}

		// set natural boundary conditions
		for idx, fc := range cell.FaceBcs {
			o.NatBcs = append(o.NatBcs, &ele.NaturalBc{fc.Cond, fc.FaceId, fc.Func, fc.Extra})

			// allocate extrapolation structures
			if fc.Cond == "ql" || fc.Cond == "seep" {
				nv := o.Cell.Shp.Nverts
				nip := len(o.IpsElem)
				o.Rhol_ex = make([]float64, nv)
				o.Drholdpl_ex = la.MatAlloc(nv, nv)
				o.Emat = la.MatAlloc(nv, nip)
				o.DoExtrap = true
				err = o.Cell.Shp.Extrapolator(o.Emat, o.IpsElem)
				if err != nil {
					chk.Panic("cannot build extrapolator matrix for p-element:\n%v", err)
				}
			}

			// additional seepage condition structures: hydrostatic flags
			if fc.Cond == "seep" {
				if len(o.Hst) == 0 {
					o.Hst = make([]bool, len(cell.FaceBcs))
				}
				if s_val, found := io.Keycode(fc.Extra, "plmax"); found {
					o.Hst[idx] = (s_val == "hst")
				}
			}
		}

		// return new element
		return &o
	})
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o *Liquid) Id() int { return o.Cell.Id }

// SetEqs sets equations
func (o *Liquid) SetEqs(eqs [][]int, mixedform_eqs []int) (err error) {
	o.Pmap = make([]int, o.Np)
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		o.Pmap[m] = eqs[m][0]
	}
	if o.HasSeep {
		for i, m := range o.SeepId2vid {
			o.Fmap[i] = eqs[m][1]
		}
	}
	return
}

// SetEleConds sets element conditions
func (o *Liquid) SetEleConds(key string, f fun.Func, extra string) (err error) {
	if key == "g" { // gravity
		o.Gfcn = f
	}
	return
}

// InterpStarVars interpolates star variables to integration points
func (o *Liquid) InterpStarVars(sol *ele.Solution) (err error) {

	// for each integration point
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}

		// interpolate starred variables
		o.PsiL[idx] = 0
		for m := 0; m < o.Cell.Shp.Nverts; m++ {
			o.PsiL[idx] += o.Cell.Shp.S[m] * sol.Psi[o.Pmap[m]]
		}
	}
	return
}

// AddToRhs adds -R to global residual vector fb
func (o *Liquid) AddToRhs(fb []float64, sol *ele.Solution) (err error) {

	// clear variables
	if o.DoExtrap {
		la.VecFill(o.Rhol_ex, 0)
	}

	// for each integration point
	O := o.LsVars
	β1 := sol.DynCfs.GetBet1()
	nverts := o.Cell.Shp.Nverts
	var coef, plt, klr, ρL, ρl, Cpl float64
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		err = o.CalcIpVars(idx, sol)
		if err != nil {
			return
		}
		coef = o.Cell.Shp.J * ip[3]
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G

		// tpm variables
		plt = β1*o.Pl - o.PsiL[idx]
		klr = o.Mdl.Cnd.Klr(o.States[idx].A_sl)
		ρL = o.States[idx].A_ρL
		err = o.Mdl.CalcLs(O, o.States[idx], o.Pl, 0, false)
		if err != nil {
			return
		}
		ρl = O.A_ρl
		Cpl = O.Cpl

		// compute ρwl. see Eq. (6) of [1]
		for i := 0; i < o.Ndim; i++ {
			o.Rhowl[i] = 0
			for j := 0; j < o.Ndim; j++ {
				o.Rhowl[i] += klr * o.Mdl.Klsat[i][j] * (ρL*o.Grav[j] - o.GradPl[j])
			}
		}

		// add negative of residual term to fb. see Eqs. (12) and (17) of [1]
		for m := 0; m < nverts; m++ {
			r := o.Pmap[m]
			fb[r] -= coef * S[m] * Cpl * plt
			for i := 0; i < o.Ndim; i++ {
				fb[r] += coef * G[m][i] * o.Rhowl[i] // += coef * div(ρl*wl)
			}
			if o.DoExtrap { // Eq. (19)
				o.Rhol_ex[m] += o.Emat[m][idx] * ρl
			}
		}
	}

	// contribution from natural boundary conditions
	if len(o.NatBcs) > 0 {
		return o.AddNatBcsToRhs(fb, sol)
	}
	return
}

// AddToKb adds element K to global Jacobian matrix Kb
func (o *Liquid) AddToKb(Kb *la.Triplet, sol *ele.Solution, firstIt bool) (err error) {

	// clear matrices
	la.MatFill(o.Kpp, 0)
	nverts := o.Cell.Shp.Nverts
	if o.DoExtrap {
		for i := 0; i < nverts; i++ {
			o.Rhol_ex[i] = 0
			for j := 0; j < nverts; j++ {
				o.Drholdpl_ex[i][j] = 0
			}
		}
	}

	// for each integration point
	O := o.LsVars
	Cl := o.Mdl.Liq.C
	β1 := sol.DynCfs.GetBet1()
	var coef, plt, klr, ρL, ρl, Cpl float64
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		err = o.CalcIpVars(idx, sol)
		if err != nil {
			return
		}
		coef = o.Cell.Shp.J * ip[3]
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G

		// tpm variables
		plt = β1*o.Pl - o.PsiL[idx]
		klr = o.Mdl.Cnd.Klr(o.States[idx].A_sl)
		ρL = o.States[idx].A_ρL
		err = o.Mdl.CalcLs(O, o.States[idx], o.Pl, 0, true)
		if err != nil {
			return
		}
		ρl = O.A_ρl
		Cpl = O.Cpl

		// Kpp := dRpl/dpl. see Eqs. (18), (A.2) and (A.3) of [1]
		for n := 0; n < nverts; n++ {
			for j := 0; j < o.Ndim; j++ {
				o.Tmp[j] = S[n]*O.Dklrdpl*(ρL*o.Grav[j]-o.GradPl[j]) + klr*(S[n]*Cl*o.Grav[j]-G[n][j])
			}
			for m := 0; m < nverts; m++ {
				o.Kpp[m][n] += coef * S[m] * S[n] * (O.DCpldpl*plt + β1*Cpl)
				for i := 0; i < o.Ndim; i++ {
					for j := 0; j < o.Ndim; j++ {
						o.Kpp[m][n] -= coef * G[m][i] * o.Mdl.Klsat[i][j] * o.Tmp[j]
					}
				}
				if o.DoExtrap { // inner summation term in Eq. (22)
					o.Drholdpl_ex[m][n] += o.Emat[m][idx] * Cpl * S[n]
				}
			}
			if o.DoExtrap { // Eq. (19)
				o.Rhol_ex[n] += o.Emat[n][idx] * ρl
			}
		}
	}

	// add to Kb
	if o.HasSeep {

		// contribution from natural boundary conditions
		err = o.AddNatBcsToJac(sol)
		if err != nil {
			return
		}

		// add to sparse matrix Kb
		for i, I := range o.Pmap {
			for j, J := range o.Pmap {
				Kb.Put(I, J, o.Kpp[i][j])
			}
			for j, J := range o.Fmap {
				Kb.Put(I, J, o.Kpf[i][j])
				Kb.Put(J, I, o.Kfp[j][i])
			}
		}
		for i, I := range o.Fmap {
			for j, J := range o.Fmap {
				Kb.Put(I, J, o.Kff[i][j])
			}
		}

	} else {

		// add to sparse matrix Kb
		for i, I := range o.Pmap {
			for j, J := range o.Pmap {
				Kb.Put(I, J, o.Kpp[i][j])
			}
		}
	}
	return
}

// Update performs (tangent) update
func (o *Liquid) Update(sol *ele.Solution) (err error) {

	// for each integration point
	var pl, Δpl float64
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, false)
		if err != nil {
			return
		}

		// compute pl and Δpl @ ip by means of interpolating from nodes
		pl, Δpl = 0, 0
		for m := 0; m < o.Cell.Shp.Nverts; m++ {
			r := o.Pmap[m]
			pl += o.Cell.Shp.S[m] * sol.Y[r]
			Δpl += o.Cell.Shp.S[m] * sol.ΔY[r]
		}

		// update state
		err = o.Mdl.Update(o.States[idx], Δpl, 0, pl, 0)
		if err != nil {
			return
		}
	}
	return
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// SetIniIvs sets initial ivs for given values in sol and ivs map
func (o *Liquid) SetIniIvs(sol *ele.Solution, ignored map[string][]float64) (err error) {

	// auxiliary
	nip := len(o.IpsElem)
	nverts := o.Cell.Shp.Nverts
	var ρL, ρG, pl, pg float64

	// allocate slices of states
	o.States = make([]*porous.State, nip)
	o.StatesBkp = make([]*porous.State, nip)
	o.StatesAux = make([]*porous.State, nip)

	// for each integration point
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G

		// interpolate pl variables
		pl = 0
		for i := 0; i < o.Ndim; i++ {
			o.GradPl[i] = 0
		}
		for m := 0; m < nverts; m++ {
			r := o.Pmap[m]
			pl += S[m] * sol.Y[r]
			for i := 0; i < o.Ndim; i++ {
				o.GradPl[i] += G[m][i] * sol.Y[r]
			}
		}

		// compute density from hydrostatic condition => enforce initial ρwl = 0
		ρL = o.Mdl.Liq.R0
		o.ComputeGrav(sol.T)
		if math.Abs(o.Grav[o.Ndim-1]) > 0 {
			ρL = o.GradPl[o.Ndim-1] / o.Grav[o.Ndim-1]
		}

		// state initialisation
		o.States[idx], err = o.Mdl.NewState(ρL, ρG, pl, pg)
		if err != nil {
			return
		}

		// backup copy
		o.StatesBkp[idx] = o.States[idx].GetCopy()
		o.StatesAux[idx] = o.States[idx].GetCopy()
	}

	// seepage face structures
	if o.HasSeep {
		o.Plmax = la.MatAlloc(len(o.NatBcs), len(o.IpsFace))
		for idx, nbc := range o.NatBcs {
			iface := nbc.IdxFace
			for jdx, ipf := range o.IpsFace {
				err = o.Cell.Shp.CalcAtFaceIp(o.X, ipf, iface)
				if err != nil {
					return
				}
				Sf := o.Cell.Shp.Sf
				switch nbc.Key {
				case "seep":
					pl = 0
					for i, m := range o.Cell.Shp.FaceLocalVerts[iface] {
						pl += Sf[i] * sol.Y[o.Pmap[m]]
					}
					o.Plmax[idx][jdx] = pl
				}
			}
		}
	}
	return
}

// BackupIvs creates copy of internal variables
func (o *Liquid) BackupIvs(aux bool) (err error) {
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

// RestoreIvs restores internal variables from copies
func (o *Liquid) RestoreIvs(aux bool) (err error) {
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
func (o *Liquid) Ureset(sol *ele.Solution) (err error) {
	return
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o *Liquid) Encode(enc utl.Encoder) (err error) {
	return enc.Encode(o.States)
}

// Decode decodes internal variables
func (o *Liquid) Decode(dec utl.Decoder) (err error) {
	err = dec.Decode(&o.States)
	if err != nil {
		return
	}
	return o.BackupIvs(false)
}

// OutIpCoords returns the coordinates of integration points
func (o *Liquid) OutIpCoords() (C [][]float64) {
	C = make([][]float64, len(o.IpsElem))
	for idx, ip := range o.IpsElem {
		C[idx] = o.Cell.Shp.IpRealCoords(o.X, ip)
	}
	return
}

// OutIpKeys returns the integration points' keys
func (o *Liquid) OutIpKeys() []string {
	return append([]string{"pl", "sl"}, LiqFlowKeys(o.Ndim)...)
}

// OutIpVals returns the integration points' values corresponding to keys
func (o *Liquid) OutIpVals(M *ele.IpsMap, sol *ele.Solution) {
	flow := LiqFlowKeys(o.Ndim)
	nip := len(o.IpsElem)
	for idx, _ := range o.IpsElem {
		err := o.CalcIpVars(idx, sol)
		if err != nil {
			return
		}
		sl := o.States[idx].A_sl
		ρL := o.States[idx].A_ρL
		klr := o.Mdl.Cnd.Klr(sl)
		M.Set("pl", idx, nip, o.Pl)
		M.Set("sl", idx, nip, sl)
		for i := 0; i < o.Ndim; i++ {
			var nwl_i float64
			for j := 0; j < o.Ndim; j++ {
				nwl_i += klr * o.Mdl.Klsat[i][j] * (o.Grav[j] - o.GradPl[j]/ρL)
			}
			M.Set(flow[i], idx, nip, nwl_i)
		}
	}
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// CalcIpVars computes current values @ integration points. idx == index of integration point
func (o *Liquid) CalcIpVars(idx int, sol *ele.Solution) (err error) {

	// interpolation functions and gradients
	err = o.Cell.Shp.CalcAtIp(o.X, o.IpsElem[idx], true)
	if err != nil {
		return
	}

	// auxiliary
	o.ComputeGrav(sol.T)

	// clear pl and its gradient @ ip
	o.Pl = 0
	for i := 0; i < o.Ndim; i++ {
		o.GradPl[i] = 0
	}

	// compute pl and its gradient @ ip by means of interpolating from nodes
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		r := o.Pmap[m]
		o.Pl += o.Cell.Shp.S[m] * sol.Y[r]
		for i := 0; i < o.Ndim; i++ {
			o.GradPl[i] += o.Cell.Shp.G[m][i] * sol.Y[r]
		}
	}
	return
}

// CalcFaceIpVars computes current values @ face integration points
func (o *Liquid) CalcFaceIpVars(fidx int, sol *ele.Solution) (ρl, pl, fl float64) {
	Sf := o.Cell.Shp.Sf
	for i, m := range o.Cell.Shp.FaceLocalVerts[fidx] {
		μ := o.Vid2seepId[m]
		ρl += Sf[i] * o.Rhol_ex[m]
		pl += Sf[i] * sol.Y[o.Pmap[m]]
		fl += Sf[i] * sol.Y[o.Fmap[μ]]
	}
	return
}

// AddNatBcsToRhs adds natural boundary conditions to rhs
func (o *Liquid) AddNatBcsToRhs(fb []float64, sol *ele.Solution) (err error) {

	// compute surface integral
	var tmp float64
	var ρl, pl, fl, plmax, g, rmp, rx, rf float64
	for idx, nbc := range o.NatBcs {

		// tmp := plmax shift or qlb
		tmp = nbc.Fcn.F(sol.T, nil)

		// loop over ips of face
		for jdx, ipf := range o.IpsFace {

			// interpolation functions and gradients @ face
			iface := nbc.IdxFace
			err = o.Cell.Shp.CalcAtFaceIp(o.X, ipf, iface)
			if err != nil {
				return
			}
			Sf := o.Cell.Shp.Sf
			Jf := la.VecNorm(o.Cell.Shp.Fnvec)
			coef := ipf[3] * Jf

			// select natural boundary condition type
			switch nbc.Key {

			// flux prescribed
			case "ql":
				ρl = 0
				for i, m := range o.Cell.Shp.FaceLocalVerts[iface] {
					ρl += Sf[i] * o.Rhol_ex[m]
				}
				for i, m := range o.Cell.Shp.FaceLocalVerts[iface] {
					fb[o.Pmap[m]] -= coef * ρl * tmp * Sf[i]
				}

			// seepage face
			case "seep":

				// variables extrapolated to face
				ρl, pl, fl = o.CalcFaceIpVars(iface, sol)
				plmax = o.Plmax[idx][jdx] - tmp
				if plmax < 0 {
					plmax = 0
				}

				// compute residuals
				g = pl - plmax // Eq. (24)
				rmp = o.Ramp(fl + o.Kap*g)
				rx = ρl * rmp // Eq. (30)
				rf = fl - rmp // Eq. (26)
				for i, m := range o.Cell.Shp.FaceLocalVerts[iface] {
					μ := o.Vid2seepId[m]
					fb[o.Pmap[m]] -= coef * Sf[i] * rx
					fb[o.Fmap[μ]] -= coef * Sf[i] * rf
				}
			}
		}
	}
	return
}

// AddNatBcsToJac adds contribution from natural boundary conditions to Jacobian
func (o *Liquid) AddNatBcsToJac(sol *ele.Solution) (err error) {

	// clear matrices
	if o.HasSeep {
		for i := 0; i < o.Np; i++ {
			for j := 0; j < o.Nf; j++ {
				o.Kpf[i][j] = 0
				o.Kfp[j][i] = 0
			}
		}
		la.MatFill(o.Kff, 0)
	}

	// compute surface integral
	nverts := o.Cell.Shp.Nverts
	var shift float64
	var ρl, pl, fl, plmax, g, rmp, rmpD float64
	var drxdpl, drxdfl, drfdpl, drfdfl float64
	for idx, nbc := range o.NatBcs {

		// plmax shift
		shift = nbc.Fcn.F(sol.T, nil)

		// loop over ips of face
		for jdx, ipf := range o.IpsFace {

			// interpolation functions and gradients @ face
			iface := nbc.IdxFace
			err = o.Cell.Shp.CalcAtFaceIp(o.X, ipf, iface)
			if err != nil {
				return
			}
			Sf := o.Cell.Shp.Sf
			Jf := la.VecNorm(o.Cell.Shp.Fnvec)
			coef := ipf[3] * Jf

			// select natural boundary condition type
			switch nbc.Key {
			case "seep":

				// variables extrapolated to face
				ρl, pl, fl = o.CalcFaceIpVars(iface, sol)
				plmax = o.Plmax[idx][jdx] - shift
				if plmax < 0 {
					plmax = 0
				}

				// compute derivatives
				g = pl - plmax // Eq. (24)
				rmp = o.Ramp(fl + o.Kap*g)
				rmpD = o.RampDeriv(fl + o.Kap*g)
				drxdpl = ρl * o.Kap * rmpD // first term in Eq. (A.4) (without Sn)
				drxdfl = ρl * rmpD         // Eq. (A.5) (without Sn)
				drfdpl = -o.Kap * rmpD     // Eq. (A.6) (corrected with κ and without Sn)
				drfdfl = 1.0 - rmpD        // Eq. (A.7) (without Sn)
				for i, m := range o.Cell.Shp.FaceLocalVerts[iface] {
					μ := o.Vid2seepId[m]
					for j, n := range o.Cell.Shp.FaceLocalVerts[iface] {
						ν := o.Vid2seepId[n]
						o.Kpp[m][n] += coef * Sf[i] * Sf[j] * drxdpl
						o.Kpf[m][ν] += coef * Sf[i] * Sf[j] * drxdfl
						o.Kfp[μ][n] += coef * Sf[i] * Sf[j] * drfdpl
						o.Kff[μ][ν] += coef * Sf[i] * Sf[j] * drfdfl
					}
					for n := 0; n < nverts; n++ { // Eqs. (18) and (22)
						for l, r := range o.Cell.Shp.FaceLocalVerts[iface] {
							o.Kpp[m][n] += coef * Sf[i] * Sf[l] * o.Drholdpl_ex[r][n] * rmp
						}
					}
				}
			}
		}
	}
	return
}

// Ramp implements the ramp function
func (o *Liquid) Ramp(x float64) float64 {
	if o.Macaulay {
		return fun.Ramp(x)
	}
	return fun.Sramp(x, o.BetRmp)
}

// RampDeriv returns the ramp function first derivative
func (o *Liquid) RampDeriv(x float64) float64 {
	if o.Macaulay {
		return fun.Heav(x)
	}
	return fun.SrampD1(x, o.BetRmp)
}

// ComputeGrav computes gravity vector @ time t
func (o *Liquid) ComputeGrav(t float64) {
	o.Grav[o.Ndim-1] = 0
	if o.Gfcn != nil {
		o.Grav[o.Ndim-1] = -o.Gfcn.F(t, nil)
	}
}
