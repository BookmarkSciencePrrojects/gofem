// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mdl/por"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// ElemPP implements an element for liquid-gas flow analyses
type ElemPP struct {

	// basic data
	Cell *inp.Cell   // the cell structure
	X    [][]float64 // matrix of nodal coordinates [ndim][nnode]
	Np   int         // number of vertices == number of pl or pg
	Ndim int         // space dimension

	// integration points
	IpsElem []shp.Ipoint // integration points of element
	IpsFace []shp.Ipoint // integration points corresponding to faces

	// material model
	Mdl *por.Model // model

	// problem variables
	Plmap []int // assembly map -- liquid pressure
	Pgmap []int // assembly map -- gas pressure

	// internal variables
	States    []*por.State
	StatesBkp []*por.State
	StatesAux []*por.State

	// gravity
	Gfcn fun.Func // gravity function

	// natural boundary conditions
	NatBcs []*NaturalBc // natural boundary conditions

	// flux boundary conditions (qb == \bar{q})
	ρl_ex     []float64   // [nverts] ρl extrapolted to nodes => if has qb (flux)
	ρg_ex     []float64   // [nverts] ρg extrapolted to nodes => if has qb (flux)
	dρldpl_ex [][]float64 // [nverts][nverts] Cpl extrapolted to nodes => if has qb (flux)
	dρgdpg_ex [][]float64 // [nverts][nverts] Cpl extrapolted to nodes => if has qb (flux)
	Emat      [][]float64 // [nverts][nips] extrapolator matrix
	DoExtrap  bool        // do extrapolation of ρl and Cpl => for use with flux and seepage conditions

	// seepage face
	Nf         int         // number of fl variables
	HasSeep    bool        // indicates if this element has seepage faces
	Vid2seepId []int       // [nverts] maps local vertex id to index in Fmap
	SeepId2vid []int       // [nf] maps seepage face variable id to local vertex id
	Flmap      []int       // [nf] map of "fl" variables (seepage face)
	Macaulay   bool        // use discrete ramp function instead of smooth ramp
	βrmp       float64     // coefficient for Sramp
	κ          float64     // κ coefficient to normalise equation for seepage face modelling
	Hst        []bool      // [nf] set hydrostatic plmax
	Plmax      [][]float64 // [nf][nipsFace] specified plmax (not corrected by multiplier)

	// local starred variables
	ψl []float64 // [nip] ψl* = β1・pl + β2・dpldt
	ψg []float64 // [nip] ψg* = β1・pg + β2・dpgdt

	// scratchpad. computed @ each ip
	g         []float64    // [ndim] gravity vector
	pl        float64      // pl: liquid pressure
	pg        float64      // pg: gas pressure
	gpl       []float64    // [ndim] ∇pl: gradient of liquid pressure
	gpg       []float64    // [ndim] ∇pg: gradient of gas pressure
	wlb       []float64    // [ndim] wlb = ρl*wl: augmented filter velocity -- liquid
	wgb       []float64    // [ndim] wgb = ρg*wg: augmented filter velocity -- gas
	dwlbdpl_n []float64    // [ndim] dwlb/dpl^n
	dwlbdpg_n []float64    // [ndim] dwlb/dpg^n
	dwgbdpl_n []float64    // [ndim] dwgb/dpl^n
	dwgbdpg_n []float64    // [ndim] dwgb/dpg^n
	Kll       [][]float64  // [np][np] dRpl/dpl consistent tangent matrix
	Klg       [][]float64  // [np][np] dRpl/dpg consistent tangent matrix
	Kgl       [][]float64  // [np][np] dRpg/dpl consistent tangent matrix
	Kgg       [][]float64  // [np][np] dRpg/dpg consistent tangent matrix
	Klf       [][]float64  // [np][nf] dRpl/dfl consistent tangent matrix
	Kfl       [][]float64  // [nf][np] dRfl/dpl consistent tangent matrix
	Kff       [][]float64  // [nf][nf] dRfl/dfl consistent tangent matrix
	res       *por.LgsVars // variable to hold results from CalcLgs
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	infogetters["pp"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *Info {

		// new info
		var info Info

		// number of nodes in element
		nverts := cell.GetNverts(edat.Lbb)

		// solution variables
		ykeys := []string{"pl", "pg"}
		info.Dofs = make([][]string, nverts)
		for m := 0; m < nverts; m++ {
			info.Dofs[m] = ykeys
		}

		// maps
		info.Y2F = map[string]string{"pl": "ql", "pg": "qg"}

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
	}

	// element allocator
	eallocators["pp"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) Elem {

		// basic data
		var o ElemPP
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
		o.ψl = make([]float64, nip)
		o.ψg = make([]float64, nip)

		// scratchpad. computed @ each ip
		o.g = make([]float64, o.Ndim)
		o.gpl = make([]float64, o.Ndim)
		o.gpg = make([]float64, o.Ndim)
		o.wlb = make([]float64, o.Ndim)
		o.wgb = make([]float64, o.Ndim)
		o.dwlbdpl_n = make([]float64, o.Ndim)
		o.dwlbdpg_n = make([]float64, o.Ndim)
		o.dwgbdpl_n = make([]float64, o.Ndim)
		o.dwgbdpg_n = make([]float64, o.Ndim)
		o.Kll = la.MatAlloc(o.Np, o.Np)
		o.Klg = la.MatAlloc(o.Np, o.Np)
		o.Kgl = la.MatAlloc(o.Np, o.Np)
		o.Kgg = la.MatAlloc(o.Np, o.Np)
		o.res = new(por.LgsVars)

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
			o.Flmap = make([]int, o.Nf)
			for μ, m := range o.SeepId2vid {
				o.Vid2seepId[m] = μ
			}

			// flags
			o.Macaulay, o.βrmp, o.κ = GetSeepFaceFlags(edat.Extra)

			// allocate coupling matrices
			o.Klf = la.MatAlloc(o.Np, o.Nf)
			o.Kfl = la.MatAlloc(o.Nf, o.Np)
			o.Kff = la.MatAlloc(o.Nf, o.Nf)
		}

		// set natural boundary conditions
		for idx, fc := range cell.FaceBcs {
			o.NatBcs = append(o.NatBcs, &NaturalBc{fc.Cond, fc.FaceId, fc.Func, fc.Extra})

			// allocate extrapolation structures
			if fc.Cond == "ql" || fc.Cond == "qg" || fc.Cond == "seep" {
				nv := o.Cell.Shp.Nverts
				nip := len(o.IpsElem)
				o.ρl_ex = make([]float64, nv)
				o.ρg_ex = make([]float64, nv)
				o.dρldpl_ex = la.MatAlloc(nv, nv)
				o.dρgdpg_ex = la.MatAlloc(nv, nv)
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
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o *ElemPP) Id() int { return o.Cell.Id }

// SetEqs sets equations
func (o *ElemPP) SetEqs(eqs [][]int, mixedform_eqs []int) (err error) {
	o.Plmap = make([]int, o.Np)
	o.Pgmap = make([]int, o.Np)
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		o.Plmap[m] = eqs[m][0]
		o.Pgmap[m] = eqs[m][1]
	}
	if o.HasSeep {
		for i, m := range o.SeepId2vid {
			o.Flmap[i] = eqs[m][2]
		}
	}
	return
}

// SetEleConds sets element conditions
func (o *ElemPP) SetEleConds(key string, f fun.Func, extra string) (err error) {
	if key == "g" { // gravity
		o.Gfcn = f
	}
	return
}

// InterpStarVars interpolates star variables to integration points
func (o *ElemPP) InterpStarVars(sol *Solution) (err error) {

	// for each integration point
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}

		// interpolate starred variables
		o.ψl[idx], o.ψg[idx] = 0, 0
		for m := 0; m < o.Cell.Shp.Nverts; m++ {
			o.ψl[idx] += o.Cell.Shp.S[m] * sol.Psi[o.Plmap[m]]
			o.ψg[idx] += o.Cell.Shp.S[m] * sol.Psi[o.Pgmap[m]]
		}
	}
	return
}

// AddToRhs adds -R to global residual vector fb
func (o *ElemPP) AddToRhs(fb []float64, sol *Solution) (err error) {

	// clear variables
	if o.DoExtrap {
		la.VecFill(o.ρl_ex, 0)
		la.VecFill(o.ρg_ex, 0)
	}

	// for each integration point
	O := o.res
	β1 := sol.DynCfs.β1
	nverts := o.Cell.Shp.Nverts
	var coef, plt, pgt, klr, kgr, ρL, ρG, ρl, ρg float64
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		err = o.ipvars(idx, sol)
		if err != nil {
			return
		}
		coef = o.Cell.Shp.J * ip[3]
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G

		// tpm variables
		plt = β1*o.pl - o.ψl[idx]
		pgt = β1*o.pg - o.ψg[idx]
		klr = o.Mdl.Cnd.Klr(o.States[idx].A_sl)
		kgr = o.Mdl.Cnd.Kgr(1.0 - o.States[idx].A_sl)
		ρL = o.States[idx].A_ρL
		ρG = o.States[idx].A_ρG
		err = o.Mdl.CalcLgs(o.res, o.States[idx], o.pl, o.pg, 0, false)
		if err != nil {
			return
		}
		ρl = o.res.A_ρl
		ρg = o.res.A_ρg

		// compute augmented filter velocities
		for i := 0; i < o.Ndim; i++ {
			o.wlb[i], o.wgb[i] = 0, 0
			for j := 0; j < o.Ndim; j++ {
				o.wlb[i] += klr * o.Mdl.Klsat[i][j] * (ρL*o.g[j] - o.gpl[j])
				o.wgb[i] += kgr * o.Mdl.Kgsat[i][j] * (ρG*o.g[j] - o.gpg[j])
			}
		}

		// add negative of residual term to fb
		for m := 0; m < nverts; m++ {
			rl := o.Plmap[m]
			rg := o.Pgmap[m]
			fb[rl] -= coef * S[m] * (O.Cpl*plt + O.Cpg*pgt)
			fb[rg] -= coef * S[m] * (O.Dpl*plt + O.Dpg*pgt)
			for i := 0; i < o.Ndim; i++ {
				fb[rl] += coef * G[m][i] * o.wlb[i] // += coef * div(ρl*wl)
				fb[rg] += coef * G[m][i] * o.wgb[i] // += coef * div(ρg*wg)
			}
			if o.DoExtrap {
				o.ρl_ex[m] += o.Emat[m][idx] * ρl
				o.ρg_ex[m] += o.Emat[m][idx] * ρg
			}
		}
	}

	// contribution from natural boundary conditions
	if len(o.NatBcs) > 0 {
		return o.add_natbcs_to_rhs(fb, sol)
	}
	return
}

// AddToKb adds element K to global Jacobian matrix Kb
func (o *ElemPP) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (err error) {

	// clear matrices
	for i := 0; i < o.Np; i++ {
		for j := 0; j < o.Np; j++ {
			o.Kll[i][j], o.Klg[i][j] = 0, 0
			o.Kgl[i][j], o.Kgg[i][j] = 0, 0
		}
	}
	nverts := o.Cell.Shp.Nverts
	if o.DoExtrap {
		for i := 0; i < nverts; i++ {
			o.ρl_ex[i], o.ρg_ex[i] = 0, 0
			for j := 0; j < nverts; j++ {
				o.dρldpl_ex[i][j] = 0
				o.dρgdpg_ex[i][j] = 0
			}
		}
	}

	// for each integration point
	O := o.res
	Cl := o.Mdl.Liq.C
	Cg := o.Mdl.Gas.C
	β1 := sol.DynCfs.β1
	var coef, plt, pgt, klr, kgr, ρL, ρG, ρl, ρg float64
	var hl_j, hg_j, dhldpl_nj, dhgdpg_nj float64
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		err = o.ipvars(idx, sol)
		if err != nil {
			return
		}
		coef = o.Cell.Shp.J * ip[3]
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G

		// tpm variables
		plt = β1*o.pl - o.ψl[idx]
		pgt = β1*o.pg - o.ψg[idx]
		klr = o.Mdl.Cnd.Klr(o.States[idx].A_sl)
		kgr = o.Mdl.Cnd.Kgr(1.0 - o.States[idx].A_sl)
		ρL = o.States[idx].A_ρL
		ρG = o.States[idx].A_ρG
		err = o.Mdl.CalcLgs(o.res, o.States[idx], o.pl, o.pg, 0, true)
		if err != nil {
			return
		}
		ρl = o.res.A_ρl
		ρg = o.res.A_ρg

		// Jacobian
		for n := 0; n < nverts; n++ {
			for i := 0; i < o.Ndim; i++ {
				o.dwlbdpl_n[i], o.dwlbdpg_n[i] = 0, 0
				o.dwgbdpl_n[i], o.dwgbdpg_n[i] = 0, 0
			}
			for j := 0; j < o.Ndim; j++ {
				hl_j = ρL*o.g[j] - o.gpl[j]
				hg_j = ρG*o.g[j] - o.gpg[j]
				dhldpl_nj = S[n]*Cl*o.g[j] - G[n][j]
				dhgdpg_nj = S[n]*Cg*o.g[j] - G[n][j]
				for i := 0; i < o.Ndim; i++ {
					o.dwlbdpl_n[i] += o.Mdl.Klsat[i][j] * (S[n]*O.Dklrdpl*hl_j + klr*dhldpl_nj)
					o.dwlbdpg_n[i] += o.Mdl.Klsat[i][j] * (S[n] * O.Dklrdpg * hl_j)
					o.dwgbdpl_n[i] += o.Mdl.Kgsat[i][j] * (S[n] * O.Dkgrdpl * hg_j)
					o.dwgbdpg_n[i] += o.Mdl.Kgsat[i][j] * (S[n]*O.Dkgrdpg*hg_j + kgr*dhgdpg_nj)
				}
			}
			for m := 0; m < nverts; m++ {
				o.Kll[m][n] += coef * S[m] * S[n] * (O.DCpldpl*plt + O.DCpgdpl*pgt + β1*O.Cpl)
				o.Klg[m][n] += coef * S[m] * S[n] * (O.DCpldpg*plt + O.DCpgdpg*pgt + β1*O.Cpg)
				o.Kgl[m][n] += coef * S[m] * S[n] * (O.DDpldpl*plt + O.DDpgdpl*pgt + β1*O.Dpl)
				o.Kgg[m][n] += coef * S[m] * S[n] * (O.DDpldpg*plt + O.DDpgdpg*pgt + β1*O.Dpg)
				for i := 0; i < o.Ndim; i++ {
					o.Kll[m][n] -= coef * G[m][i] * o.dwlbdpl_n[i]
					o.Klg[m][n] -= coef * G[m][i] * o.dwlbdpg_n[i]
					o.Kgl[m][n] -= coef * G[m][i] * o.dwgbdpl_n[i]
					o.Kgg[m][n] -= coef * G[m][i] * o.dwgbdpg_n[i]
				}
				if o.DoExtrap {
					o.dρldpl_ex[m][n] += o.Emat[m][idx] * O.Cpl * S[n]
					o.dρgdpg_ex[m][n] += o.Emat[m][idx] * O.Dpg * S[n]
				}
			}
			if o.DoExtrap {
				o.ρl_ex[n] += o.Emat[n][idx] * ρl
				o.ρg_ex[n] += o.Emat[n][idx] * ρg
			}
		}
	}

	// contribution from natural boundary conditions
	if o.HasSeep {
		err = o.add_natbcs_to_jac(sol)
		if err != nil {
			return
		}
	}

	// assemble K matrices into Kb
	o.assembleKs(Kb)
	return
}

// assembleKs assemble K matrices into Kb
func (o *ElemPP) assembleKs(Kb *la.Triplet) {
	if o.HasSeep {
		for i, I := range o.Plmap {
			for j, J := range o.Plmap {
				Kb.Put(I, J, o.Kll[i][j])
			}
			for j, J := range o.Pgmap {
				Kb.Put(I, J, o.Klg[i][j])
				Kb.Put(J, I, o.Kgl[j][i])
			}
			for j, J := range o.Flmap {
				Kb.Put(I, J, o.Klf[i][j])
				Kb.Put(J, I, o.Kfl[j][i])
			}
		}
		for i, I := range o.Pgmap {
			for j, J := range o.Pgmap {
				Kb.Put(I, J, o.Kgg[i][j])
			}
		}
		for i, I := range o.Flmap {
			for j, J := range o.Flmap {
				Kb.Put(I, J, o.Kff[i][j])
			}
		}
	} else {
		for i, I := range o.Plmap {
			for j, J := range o.Plmap {
				Kb.Put(I, J, o.Kll[i][j])
			}
			for j, J := range o.Pgmap {
				Kb.Put(I, J, o.Klg[i][j])
				Kb.Put(J, I, o.Kgl[j][i])
			}
		}
		for i, I := range o.Pgmap {
			for j, J := range o.Pgmap {
				Kb.Put(I, J, o.Kgg[i][j])
			}
		}
	}
}

// Update performs (tangent) update
func (o *ElemPP) Update(sol *Solution) (err error) {

	// for each integration point
	var pl, pg, Δpl, Δpg float64
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, false)
		if err != nil {
			return
		}

		// compute pl, pg, Δpl and Δpg @ ip by means of interpolating from nodes
		pl, pg, Δpl, Δpg = 0, 0, 0, 0
		for m := 0; m < o.Cell.Shp.Nverts; m++ {
			rl := o.Plmap[m]
			rg := o.Pgmap[m]
			pl += o.Cell.Shp.S[m] * sol.Y[rl]
			pg += o.Cell.Shp.S[m] * sol.Y[rg]
			Δpl += o.Cell.Shp.S[m] * sol.ΔY[rl]
			Δpg += o.Cell.Shp.S[m] * sol.ΔY[rg]
		}

		// update state
		err = o.Mdl.Update(o.States[idx], Δpl, Δpg, pl, pg)
		if err != nil {
			return
		}
	}
	return
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// SetIniIvs sets initial ivs for given values in sol and ivs map
//  Note: this function assumes a hydrostatic fully saturated initial condition, thus:
//        sl=1, krl=1, wlb=0  sg=0, krg=krgMin, wgb=0
func (o *ElemPP) SetIniIvs(sol *Solution, ignored map[string][]float64) (err error) {

	// auxiliary
	nip := len(o.IpsElem)
	nverts := o.Cell.Shp.Nverts
	var ρL, ρG, pl, pg float64

	// allocate slices of states
	o.States = make([]*por.State, nip)
	o.StatesBkp = make([]*por.State, nip)
	o.StatesAux = make([]*por.State, nip)

	// for each integration point
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G

		// interpolate pl and pg variables
		pl, pg = 0, 0
		for i := 0; i < o.Ndim; i++ {
			o.gpl[i], o.gpg[i] = 0, 0
		}
		for m := 0; m < nverts; m++ {
			rl := o.Plmap[m]
			rg := o.Pgmap[m]
			pl += S[m] * sol.Y[rl]
			pg += S[m] * sol.Y[rg]
			for i := 0; i < o.Ndim; i++ {
				o.gpl[i] += G[m][i] * sol.Y[rl]
				o.gpg[i] += G[m][i] * sol.Y[rg]
			}
		}

		// compute density from hydrostatic condition => enforce initial ρwl = 0
		ρL = o.Mdl.Liq.R0
		ρG = o.Mdl.Gas.R0
		o.compute_gvec(sol.T)
		if math.Abs(o.g[o.Ndim-1]) > 0 {
			ρL = o.gpl[o.Ndim-1] / o.g[o.Ndim-1]
			ρG = o.gpg[o.Ndim-1] / o.g[o.Ndim-1]
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
						pl += Sf[i] * sol.Y[o.Plmap[m]]
					}
					o.Plmax[idx][jdx] = pl
				}
			}
		}
	}
	return
}

// BackupIvs creates copy of internal variables
func (o *ElemPP) BackupIvs(aux bool) (err error) {
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
func (o *ElemPP) RestoreIvs(aux bool) (err error) {
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
func (o *ElemPP) Ureset(sol *Solution) (err error) {
	return
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o *ElemPP) Encode(enc Encoder) (err error) {
	return enc.Encode(o.States)
}

// Decode decodes internal variables
func (o *ElemPP) Decode(dec Decoder) (err error) {
	err = dec.Decode(&o.States)
	if err != nil {
		return
	}
	return o.BackupIvs(false)
}

// OutIpCoords returns the coordinates of integration points
func (o *ElemPP) OutIpCoords() (C [][]float64) {
	C = make([][]float64, len(o.IpsElem))
	for idx, ip := range o.IpsElem {
		C[idx] = o.Cell.Shp.IpRealCoords(o.X, ip)
	}
	return
}

// OutIpKeys returns the integration points' keys
func (o *ElemPP) OutIpKeys() []string {
	keys := append([]string{"pl", "pg", "sl"}, LiqFlowKeys(o.Ndim)...)
	return append(keys, GasFlowKeys(o.Ndim)...)
}

// OutIpVals returns the integration points' values corresponding to keys
func (o *ElemPP) OutIpVals(M *IpsMap, sol *Solution) {
	flowL := LiqFlowKeys(o.Ndim)
	flowG := GasFlowKeys(o.Ndim)
	nip := len(o.IpsElem)
	for idx, _ := range o.IpsElem {
		err := o.ipvars(idx, sol)
		if err != nil {
			return
		}
		sl := o.States[idx].A_sl
		sg := 1.0 - sl
		ρL := o.States[idx].A_ρL
		ρG := o.States[idx].A_ρG
		klr := o.Mdl.Cnd.Klr(sl)
		kgr := o.Mdl.Cnd.Klr(sg)
		M.Set("pl", idx, nip, o.pl)
		M.Set("pg", idx, nip, o.pg)
		M.Set("sl", idx, nip, sl)
		for i := 0; i < o.Ndim; i++ {
			var nwl_i, nwg_i float64
			for j := 0; j < o.Ndim; j++ {
				nwl_i += klr * o.Mdl.Klsat[i][j] * (o.g[j] - o.gpl[j]/ρL)
				nwg_i += kgr * o.Mdl.Kgsat[i][j] * (o.g[j] - o.gpg[j]/ρG)
			}
			M.Set(flowL[i], idx, nip, nwl_i)
			M.Set(flowG[i], idx, nip, nwg_i)
		}
	}
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// ipvars computes current values @ integration points. idx == index of integration point
func (o *ElemPP) ipvars(idx int, sol *Solution) (err error) {

	// interpolation functions and gradients
	err = o.Cell.Shp.CalcAtIp(o.X, o.IpsElem[idx], true)
	if err != nil {
		return
	}

	// auxiliary
	o.compute_gvec(sol.T)

	// clear pl, pg and gradients @ ip
	o.pl, o.pg = 0, 0
	for i := 0; i < o.Ndim; i++ {
		o.gpl[i], o.gpg[i] = 0, 0
	}

	// compute pl, pg and gradients @ ip by means of interpolating from nodes
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		rl := o.Plmap[m]
		rg := o.Pgmap[m]
		o.pl += o.Cell.Shp.S[m] * sol.Y[rl]
		o.pg += o.Cell.Shp.S[m] * sol.Y[rg]
		for i := 0; i < o.Ndim; i++ {
			o.gpl[i] += o.Cell.Shp.G[m][i] * sol.Y[rl]
			o.gpg[i] += o.Cell.Shp.G[m][i] * sol.Y[rg]
		}
	}
	return
}

// fipvars computes current values @ face integration points
func (o *ElemPP) fipvars(fidx int, sol *Solution) (ρl, pl, fl float64) {
	Sf := o.Cell.Shp.Sf
	for i, m := range o.Cell.Shp.FaceLocalVerts[fidx] {
		μ := o.Vid2seepId[m]
		ρl += Sf[i] * o.ρl_ex[m]
		pl += Sf[i] * sol.Y[o.Plmap[m]]
		fl += Sf[i] * sol.Y[o.Flmap[μ]]
	}
	return
}

// add_natbcs_to_rhs adds natural boundary conditions to rhs
func (o *ElemPP) add_natbcs_to_rhs(fb []float64, sol *Solution) (err error) {

	// compute surface integral
	var tmp float64
	var ρl, ρg, pl, fl, plmax, g, rmp, rx, rf float64
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

			// liquid flux prescribed
			case "ql":
				ρl = 0
				for i, m := range o.Cell.Shp.FaceLocalVerts[iface] {
					ρl += Sf[i] * o.ρl_ex[m]
				}
				for i, m := range o.Cell.Shp.FaceLocalVerts[iface] {
					fb[o.Plmap[m]] -= coef * ρl * tmp * Sf[i]
				}

			// gas flux prescribed
			case "qg":
				ρg = 0
				for i, m := range o.Cell.Shp.FaceLocalVerts[iface] {
					ρg += Sf[i] * o.ρg_ex[m]
				}
				for i, m := range o.Cell.Shp.FaceLocalVerts[iface] {
					fb[o.Pgmap[m]] -= coef * ρg * tmp * Sf[i]
				}

			// seepage face
			case "seep":

				// variables extrapolated to face
				ρl, pl, fl = o.fipvars(iface, sol)
				plmax = o.Plmax[idx][jdx] - tmp
				if plmax < 0 {
					plmax = 0
				}

				// compute residuals
				g = pl - plmax // Eq. (24)
				rmp = o.ramp(fl + o.κ*g)
				rx = ρl * rmp // Eq. (30)
				rf = fl - rmp // Eq. (26)
				for i, m := range o.Cell.Shp.FaceLocalVerts[iface] {
					μ := o.Vid2seepId[m]
					fb[o.Plmap[m]] -= coef * Sf[i] * rx
					fb[o.Flmap[μ]] -= coef * Sf[i] * rf
				}
			}
		}
	}
	return
}

// add_natbcs_to_jac adds contribution from natural boundary conditions to Jacobian
func (o *ElemPP) add_natbcs_to_jac(sol *Solution) (err error) {

	// clear matrices
	if o.HasSeep {
		for i := 0; i < o.Np; i++ {
			for j := 0; j < o.Nf; j++ {
				o.Klf[i][j] = 0
				o.Kfl[j][i] = 0
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
				ρl, pl, fl = o.fipvars(iface, sol)
				plmax = o.Plmax[idx][jdx] - shift
				if plmax < 0 {
					plmax = 0
				}

				// compute derivatives
				g = pl - plmax // Eq. (24)
				rmp = o.ramp(fl + o.κ*g)
				rmpD = o.rampD1(fl + o.κ*g)
				drxdpl = ρl * o.κ * rmpD // first term in Eq. (A.4) (without Sn)
				drxdfl = ρl * rmpD       // Eq. (A.5) (without Sn)
				drfdpl = -o.κ * rmpD     // Eq. (A.6) (corrected with κ and without Sn)
				drfdfl = 1.0 - rmpD      // Eq. (A.7) (without Sn)
				for i, m := range o.Cell.Shp.FaceLocalVerts[iface] {
					μ := o.Vid2seepId[m]
					for j, n := range o.Cell.Shp.FaceLocalVerts[iface] {
						ν := o.Vid2seepId[n]
						o.Kll[m][n] += coef * Sf[i] * Sf[j] * drxdpl
						o.Klf[m][ν] += coef * Sf[i] * Sf[j] * drxdfl
						o.Kfl[μ][n] += coef * Sf[i] * Sf[j] * drfdpl
						o.Kff[μ][ν] += coef * Sf[i] * Sf[j] * drfdfl
					}
					for n := 0; n < nverts; n++ { // Eqs. (18) and (22)
						for l, r := range o.Cell.Shp.FaceLocalVerts[iface] {
							o.Kll[m][n] += coef * Sf[i] * Sf[l] * o.dρldpl_ex[r][n] * rmp
						}
					}
				}
			}
		}
	}
	return
}

// ramp implements the ramp function
func (o *ElemPP) ramp(x float64) float64 {
	if o.Macaulay {
		return fun.Ramp(x)
	}
	return fun.Sramp(x, o.βrmp)
}

// rampderiv returns the ramp function first derivative
func (o *ElemPP) rampD1(x float64) float64 {
	if o.Macaulay {
		return fun.Heav(x)
	}
	return fun.SrampD1(x, o.βrmp)
}

// compute_gvec computes gravity vector @ time t
func (o *ElemPP) compute_gvec(t float64) {
	o.g[o.Ndim-1] = 0
	if o.Gfcn != nil {
		o.g[o.Ndim-1] = -o.Gfcn.F(t, nil)
	}
}
