// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package seepage

import (
	"math"

	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mdl/por"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// LiquidGas implements an element for liquid-gas flow analyses
type LiquidGas struct {

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
	NatBcs []*ele.NaturalBc // natural boundary conditions

	// flux boundary conditions (qb == \bar{q})
	Rhol_ex     []float64   // [nverts] ρl extrapolted to nodes => if has qb (flux)
	Rhog_ex     []float64   // [nverts] ρg extrapolted to nodes => if has qb (flux)
	Drholdpl_ex [][]float64 // [nverts][nverts] Cpl extrapolted to nodes => if has qb (flux)
	Drhogdpg_ex [][]float64 // [nverts][nverts] Cpl extrapolted to nodes => if has qb (flux)
	Emat        [][]float64 // [nverts][nips] extrapolator matrix
	DoExtrap    bool        // do extrapolation of ρl and Cpl => for use with flux and seepage conditions

	// seepage face
	Nf         int         // number of fl variables
	HasSeep    bool        // indicates if this element has seepage faces
	Vid2seepId []int       // [nverts] maps local vertex id to index in Fmap
	SeepId2vid []int       // [nf] maps seepage face variable id to local vertex id
	Flmap      []int       // [nf] map of "fl" variables (seepage face)
	Macaulay   bool        // use discrete ramp function instead of smooth ramp
	BetRmp     float64     // coefficient for Sramp
	Kap        float64     // Kap coefficient to normalise equation for seepage face modelling
	Hst        []bool      // [nf] set hydrostatic plmax
	Plmax      [][]float64 // [nf][nipsFace] specified plmax (not corrected by multiplier)

	// local starred variables
	PsiL []float64 // [nip] ψl* = β1・pl + β2・dpldt
	PsiG []float64 // [nip] ψg* = β1・pg + β2・dpgdt

	// scratchpad. computed @ each ip
	Grav      []float64    // [ndim] gravity vector
	Pl        float64      // pl: liquid pressure
	Pg        float64      // pg: gas pressure
	GradPl    []float64    // [ndim] ∇pl: gradient of liquid pressure
	GradPg    []float64    // [ndim] ∇pg: gradient of gas pressure
	Wlb       []float64    // [ndim] wlb = ρl*wl: augmented filter velocity -- liquid
	Wgb       []float64    // [ndim] wgb = ρg*wg: augmented filter velocity -- gas
	Dwlbdpl_n []float64    // [ndim] dwlb/dpl^n
	Dwlbdpg_n []float64    // [ndim] dwlb/dpg^n
	Dwgbdpl_n []float64    // [ndim] dwgb/dpl^n
	Dwgbdpg_n []float64    // [ndim] dwgb/dpg^n
	Kll       [][]float64  // [np][np] dRpl/dpl consistent tangent matrix
	Klg       [][]float64  // [np][np] dRpl/dpg consistent tangent matrix
	Kgl       [][]float64  // [np][np] dRpg/dpl consistent tangent matrix
	Kgg       [][]float64  // [np][np] dRpg/dpg consistent tangent matrix
	Klf       [][]float64  // [np][nf] dRpl/dfl consistent tangent matrix
	Kfl       [][]float64  // [nf][np] dRfl/dpl consistent tangent matrix
	Kff       [][]float64  // [nf][nf] dRfl/dfl consistent tangent matrix
	LgsVars   *por.LgsVars // variable to hold results from CalcLgs
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	ele.SetInfoFunc("liquid-gas", func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *ele.Info {

		// new info
		var info ele.Info

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
	})

	// element allocator
	ele.SetAllocator("liquid-gas", func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) ele.Element {

		// basic data
		var o LiquidGas
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
		o.PsiG = make([]float64, nip)

		// scratchpad. computed @ each ip
		o.Grav = make([]float64, o.Ndim)
		o.GradPl = make([]float64, o.Ndim)
		o.GradPg = make([]float64, o.Ndim)
		o.Wlb = make([]float64, o.Ndim)
		o.Wgb = make([]float64, o.Ndim)
		o.Dwlbdpl_n = make([]float64, o.Ndim)
		o.Dwlbdpg_n = make([]float64, o.Ndim)
		o.Dwgbdpl_n = make([]float64, o.Ndim)
		o.Dwgbdpg_n = make([]float64, o.Ndim)
		o.Kll = la.MatAlloc(o.Np, o.Np)
		o.Klg = la.MatAlloc(o.Np, o.Np)
		o.Kgl = la.MatAlloc(o.Np, o.Np)
		o.Kgg = la.MatAlloc(o.Np, o.Np)
		o.LgsVars = new(por.LgsVars)

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
			o.Macaulay, o.BetRmp, o.Kap = GetSeepFaceFlags(edat.Extra)

			// allocate coupling matrices
			o.Klf = la.MatAlloc(o.Np, o.Nf)
			o.Kfl = la.MatAlloc(o.Nf, o.Np)
			o.Kff = la.MatAlloc(o.Nf, o.Nf)
		}

		// set natural boundary conditions
		for idx, fc := range cell.FaceBcs {
			o.NatBcs = append(o.NatBcs, &ele.NaturalBc{fc.Cond, fc.FaceId, fc.Func, fc.Extra})

			// allocate extrapolation structures
			if fc.Cond == "ql" || fc.Cond == "qg" || fc.Cond == "seep" {
				nv := o.Cell.Shp.Nverts
				nip := len(o.IpsElem)
				o.Rhol_ex = make([]float64, nv)
				o.Rhog_ex = make([]float64, nv)
				o.Drholdpl_ex = la.MatAlloc(nv, nv)
				o.Drhogdpg_ex = la.MatAlloc(nv, nv)
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
func (o *LiquidGas) Id() int { return o.Cell.Id }

// SetEqs sets equations
func (o *LiquidGas) SetEqs(eqs [][]int, mixedform_eqs []int) (err error) {
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
func (o *LiquidGas) SetEleConds(key string, f fun.Func, extra string) (err error) {
	if key == "g" { // gravity
		o.Gfcn = f
	}
	return
}

// InterpStarVars interpolates star variables to integration points
func (o *LiquidGas) InterpStarVars(sol *ele.Solution) (err error) {

	// for each integration point
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}

		// interpolate starred variables
		o.PsiL[idx], o.PsiG[idx] = 0, 0
		for m := 0; m < o.Cell.Shp.Nverts; m++ {
			o.PsiL[idx] += o.Cell.Shp.S[m] * sol.Psi[o.Plmap[m]]
			o.PsiG[idx] += o.Cell.Shp.S[m] * sol.Psi[o.Pgmap[m]]
		}
	}
	return
}

// AddToRhs adds -R to global residual vector fb
func (o *LiquidGas) AddToRhs(fb []float64, sol *ele.Solution) (err error) {

	// clear variables
	if o.DoExtrap {
		la.VecFill(o.Rhol_ex, 0)
		la.VecFill(o.Rhog_ex, 0)
	}

	// for each integration point
	O := o.LgsVars
	β1 := sol.DynCfs.GetBet1()
	nverts := o.Cell.Shp.Nverts
	var coef, plt, pgt, klr, kgr, ρL, ρG, ρl, ρg float64
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
		pgt = β1*o.Pg - o.PsiG[idx]
		klr = o.Mdl.Cnd.Klr(o.States[idx].A_sl)
		kgr = o.Mdl.Cnd.Kgr(1.0 - o.States[idx].A_sl)
		ρL = o.States[idx].A_ρL
		ρG = o.States[idx].A_ρG
		err = o.Mdl.CalcLgs(O, o.States[idx], o.Pl, o.Pg, 0, false)
		if err != nil {
			return
		}
		ρl = O.A_ρl
		ρg = O.A_ρg

		// compute augmented filter velocities
		for i := 0; i < o.Ndim; i++ {
			o.Wlb[i], o.Wgb[i] = 0, 0
			for j := 0; j < o.Ndim; j++ {
				o.Wlb[i] += klr * o.Mdl.Klsat[i][j] * (ρL*o.Grav[j] - o.GradPl[j])
				o.Wgb[i] += kgr * o.Mdl.Kgsat[i][j] * (ρG*o.Grav[j] - o.GradPg[j])
			}
		}

		// add negative of residual term to fb
		for m := 0; m < nverts; m++ {
			rl := o.Plmap[m]
			rg := o.Pgmap[m]
			fb[rl] -= coef * S[m] * (O.Cpl*plt + O.Cpg*pgt)
			fb[rg] -= coef * S[m] * (O.Dpl*plt + O.Dpg*pgt)
			for i := 0; i < o.Ndim; i++ {
				fb[rl] += coef * G[m][i] * o.Wlb[i] // += coef * div(ρl*wl)
				fb[rg] += coef * G[m][i] * o.Wgb[i] // += coef * div(ρg*wg)
			}
			if o.DoExtrap {
				o.Rhol_ex[m] += o.Emat[m][idx] * ρl
				o.Rhog_ex[m] += o.Emat[m][idx] * ρg
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
func (o *LiquidGas) AddToKb(Kb *la.Triplet, sol *ele.Solution, firstIt bool) (err error) {

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
			o.Rhol_ex[i], o.Rhog_ex[i] = 0, 0
			for j := 0; j < nverts; j++ {
				o.Drholdpl_ex[i][j] = 0
				o.Drhogdpg_ex[i][j] = 0
			}
		}
	}

	// for each integration point
	O := o.LgsVars
	Cl := o.Mdl.Liq.C
	Cg := o.Mdl.Gas.C
	β1 := sol.DynCfs.GetBet1()
	var coef, plt, pgt, klr, kgr, ρL, ρG, ρl, ρg float64
	var hl_j, hg_j, dhldpl_nj, dhgdpg_nj float64
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
		pgt = β1*o.Pg - o.PsiG[idx]
		klr = o.Mdl.Cnd.Klr(o.States[idx].A_sl)
		kgr = o.Mdl.Cnd.Kgr(1.0 - o.States[idx].A_sl)
		ρL = o.States[idx].A_ρL
		ρG = o.States[idx].A_ρG
		err = o.Mdl.CalcLgs(O, o.States[idx], o.Pl, o.Pg, 0, true)
		if err != nil {
			return
		}
		ρl = O.A_ρl
		ρg = O.A_ρg

		// Jacobian
		for n := 0; n < nverts; n++ {
			for i := 0; i < o.Ndim; i++ {
				o.Dwlbdpl_n[i], o.Dwlbdpg_n[i] = 0, 0
				o.Dwgbdpl_n[i], o.Dwgbdpg_n[i] = 0, 0
			}
			for j := 0; j < o.Ndim; j++ {
				hl_j = ρL*o.Grav[j] - o.GradPl[j]
				hg_j = ρG*o.Grav[j] - o.GradPg[j]
				dhldpl_nj = S[n]*Cl*o.Grav[j] - G[n][j]
				dhgdpg_nj = S[n]*Cg*o.Grav[j] - G[n][j]
				for i := 0; i < o.Ndim; i++ {
					o.Dwlbdpl_n[i] += o.Mdl.Klsat[i][j] * (S[n]*O.Dklrdpl*hl_j + klr*dhldpl_nj)
					o.Dwlbdpg_n[i] += o.Mdl.Klsat[i][j] * (S[n] * O.Dklrdpg * hl_j)
					o.Dwgbdpl_n[i] += o.Mdl.Kgsat[i][j] * (S[n] * O.Dkgrdpl * hg_j)
					o.Dwgbdpg_n[i] += o.Mdl.Kgsat[i][j] * (S[n]*O.Dkgrdpg*hg_j + kgr*dhgdpg_nj)
				}
			}
			for m := 0; m < nverts; m++ {
				o.Kll[m][n] += coef * S[m] * S[n] * (O.DCpldpl*plt + O.DCpgdpl*pgt + β1*O.Cpl)
				o.Klg[m][n] += coef * S[m] * S[n] * (O.DCpldpg*plt + O.DCpgdpg*pgt + β1*O.Cpg)
				o.Kgl[m][n] += coef * S[m] * S[n] * (O.DDpldpl*plt + O.DDpgdpl*pgt + β1*O.Dpl)
				o.Kgg[m][n] += coef * S[m] * S[n] * (O.DDpldpg*plt + O.DDpgdpg*pgt + β1*O.Dpg)
				for i := 0; i < o.Ndim; i++ {
					o.Kll[m][n] -= coef * G[m][i] * o.Dwlbdpl_n[i]
					o.Klg[m][n] -= coef * G[m][i] * o.Dwlbdpg_n[i]
					o.Kgl[m][n] -= coef * G[m][i] * o.Dwgbdpl_n[i]
					o.Kgg[m][n] -= coef * G[m][i] * o.Dwgbdpg_n[i]
				}
				if o.DoExtrap {
					o.Drholdpl_ex[m][n] += o.Emat[m][idx] * O.Cpl * S[n]
					o.Drhogdpg_ex[m][n] += o.Emat[m][idx] * O.Dpg * S[n]
				}
			}
			if o.DoExtrap {
				o.Rhol_ex[n] += o.Emat[n][idx] * ρl
				o.Rhog_ex[n] += o.Emat[n][idx] * ρg
			}
		}
	}

	// contribution from natural boundary conditions
	if o.HasSeep {
		err = o.AddNatBcsToJac(sol)
		if err != nil {
			return
		}
	}

	// assemble K matrices into Kb
	o.AssembleKs(Kb)
	return
}

// AssembleKs assemble K matrices into Kb
func (o *LiquidGas) AssembleKs(Kb *la.Triplet) {
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
func (o *LiquidGas) Update(sol *ele.Solution) (err error) {

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
func (o *LiquidGas) SetIniIvs(sol *ele.Solution, ignored map[string][]float64) (err error) {

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
			o.GradPl[i], o.GradPg[i] = 0, 0
		}
		for m := 0; m < nverts; m++ {
			rl := o.Plmap[m]
			rg := o.Pgmap[m]
			pl += S[m] * sol.Y[rl]
			pg += S[m] * sol.Y[rg]
			for i := 0; i < o.Ndim; i++ {
				o.GradPl[i] += G[m][i] * sol.Y[rl]
				o.GradPg[i] += G[m][i] * sol.Y[rg]
			}
		}

		// compute density from hydrostatic condition => enforce initial ρwl = 0
		ρL = o.Mdl.Liq.R0
		ρG = o.Mdl.Gas.R0
		o.ComputeGrav(sol.T)
		if math.Abs(o.Grav[o.Ndim-1]) > 0 {
			ρL = o.GradPl[o.Ndim-1] / o.Grav[o.Ndim-1]
			ρG = o.GradPg[o.Ndim-1] / o.Grav[o.Ndim-1]
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
func (o *LiquidGas) BackupIvs(aux bool) (err error) {
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
func (o *LiquidGas) RestoreIvs(aux bool) (err error) {
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
func (o *LiquidGas) Ureset(sol *ele.Solution) (err error) {
	return
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o *LiquidGas) Encode(enc utl.Encoder) (err error) {
	return enc.Encode(o.States)
}

// Decode decodes internal variables
func (o *LiquidGas) Decode(dec utl.Decoder) (err error) {
	err = dec.Decode(&o.States)
	if err != nil {
		return
	}
	return o.BackupIvs(false)
}

// OutIpCoords returns the coordinates of integration points
func (o *LiquidGas) OutIpCoords() (C [][]float64) {
	C = make([][]float64, len(o.IpsElem))
	for idx, ip := range o.IpsElem {
		C[idx] = o.Cell.Shp.IpRealCoords(o.X, ip)
	}
	return
}

// OutIpKeys returns the integration points' keys
func (o *LiquidGas) OutIpKeys() []string {
	keys := append([]string{"pl", "pg", "sl"}, LiqFlowKeys(o.Ndim)...)
	return append(keys, GasFlowKeys(o.Ndim)...)
}

// OutIpVals returns the integration points' values corresponding to keys
func (o *LiquidGas) OutIpVals(M *ele.IpsMap, sol *ele.Solution) {
	flowL := LiqFlowKeys(o.Ndim)
	flowG := GasFlowKeys(o.Ndim)
	nip := len(o.IpsElem)
	for idx, _ := range o.IpsElem {
		err := o.CalcIpVars(idx, sol)
		if err != nil {
			return
		}
		sl := o.States[idx].A_sl
		sg := 1.0 - sl
		ρL := o.States[idx].A_ρL
		ρG := o.States[idx].A_ρG
		klr := o.Mdl.Cnd.Klr(sl)
		kgr := o.Mdl.Cnd.Klr(sg)
		M.Set("pl", idx, nip, o.Pl)
		M.Set("pg", idx, nip, o.Pg)
		M.Set("sl", idx, nip, sl)
		for i := 0; i < o.Ndim; i++ {
			var nwl_i, nwg_i float64
			for j := 0; j < o.Ndim; j++ {
				nwl_i += klr * o.Mdl.Klsat[i][j] * (o.Grav[j] - o.GradPl[j]/ρL)
				nwg_i += kgr * o.Mdl.Kgsat[i][j] * (o.Grav[j] - o.GradPg[j]/ρG)
			}
			M.Set(flowL[i], idx, nip, nwl_i)
			M.Set(flowG[i], idx, nip, nwg_i)
		}
	}
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// CalcIpVars computes current values @ integration points. idx == index of integration point
func (o *LiquidGas) CalcIpVars(idx int, sol *ele.Solution) (err error) {

	// interpolation functions and gradients
	err = o.Cell.Shp.CalcAtIp(o.X, o.IpsElem[idx], true)
	if err != nil {
		return
	}

	// auxiliary
	o.ComputeGrav(sol.T)

	// clear pl, pg and gradients @ ip
	o.Pl, o.Pg = 0, 0
	for i := 0; i < o.Ndim; i++ {
		o.GradPl[i], o.GradPg[i] = 0, 0
	}

	// compute pl, pg and gradients @ ip by means of interpolating from nodes
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		rl := o.Plmap[m]
		rg := o.Pgmap[m]
		o.Pl += o.Cell.Shp.S[m] * sol.Y[rl]
		o.Pg += o.Cell.Shp.S[m] * sol.Y[rg]
		for i := 0; i < o.Ndim; i++ {
			o.GradPl[i] += o.Cell.Shp.G[m][i] * sol.Y[rl]
			o.GradPg[i] += o.Cell.Shp.G[m][i] * sol.Y[rg]
		}
	}
	return
}

// CalcFaceIpVars computes current values @ face integration points
func (o *LiquidGas) CalcFaceIpVars(fidx int, sol *ele.Solution) (ρl, pl, fl float64) {
	Sf := o.Cell.Shp.Sf
	for i, m := range o.Cell.Shp.FaceLocalVerts[fidx] {
		μ := o.Vid2seepId[m]
		ρl += Sf[i] * o.Rhol_ex[m]
		pl += Sf[i] * sol.Y[o.Plmap[m]]
		fl += Sf[i] * sol.Y[o.Flmap[μ]]
	}
	return
}

// AddNatBcsToRhs adds natural boundary conditions to rhs
func (o *LiquidGas) AddNatBcsToRhs(fb []float64, sol *ele.Solution) (err error) {

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
					ρl += Sf[i] * o.Rhol_ex[m]
				}
				for i, m := range o.Cell.Shp.FaceLocalVerts[iface] {
					fb[o.Plmap[m]] -= coef * ρl * tmp * Sf[i]
				}

			// gas flux prescribed
			case "qg":
				ρg = 0
				for i, m := range o.Cell.Shp.FaceLocalVerts[iface] {
					ρg += Sf[i] * o.Rhog_ex[m]
				}
				for i, m := range o.Cell.Shp.FaceLocalVerts[iface] {
					fb[o.Pgmap[m]] -= coef * ρg * tmp * Sf[i]
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
					fb[o.Plmap[m]] -= coef * Sf[i] * rx
					fb[o.Flmap[μ]] -= coef * Sf[i] * rf
				}
			}
		}
	}
	return
}

// AddNatBcsToJac adds contribution from natural boundary conditions to Jacobian
func (o *LiquidGas) AddNatBcsToJac(sol *ele.Solution) (err error) {

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
				drfdpl = -o.Kap * rmpD     // Eq. (A.6) (corrected with Kap and without Sn)
				drfdfl = 1.0 - rmpD        // Eq. (A.7) (without Sn)
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
							o.Kll[m][n] += coef * Sf[i] * Sf[l] * o.Drholdpl_ex[r][n] * rmp
						}
					}
				}
			}
		}
	}
	return
}

// Ramp implements the ramp function
func (o *LiquidGas) Ramp(x float64) float64 {
	if o.Macaulay {
		return fun.Ramp(x)
	}
	return fun.Sramp(x, o.BetRmp)
}

// RampDeriv returns the ramp function first derivative
func (o *LiquidGas) RampDeriv(x float64) float64 {
	if o.Macaulay {
		return fun.Heav(x)
	}
	return fun.SrampD1(x, o.BetRmp)
}

// ComputeGrav computes gravity vector @ time t
func (o *LiquidGas) ComputeGrav(t float64) {
	o.Grav[o.Ndim-1] = 0
	if o.Gfcn != nil {
		o.Grav[o.Ndim-1] = -o.Gfcn.F(t, nil)
	}
}
