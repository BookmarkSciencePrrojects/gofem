// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package porous

import (
	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/ele/seepage"
	"github.com/cpmech/gofem/ele/solid"
	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun/dbf"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/tsr"
	"github.com/cpmech/gosl/utl"
)

// SolidLiquidGas implements the u-pl-pg formulation
type SolidLiquidGas struct {

	// auxiliary
	Sim     *inp.Simulation // simulation
	Cell    *inp.Cell       // cell
	LbbCell *inp.Cell       // if LBB==false, same as Cell; otherwise LbbCell is a new cell with less vertices
	Edat    *inp.ElemData   // element data; stored in allocator to be used in Connect
	Ndim    int             // space dimension

	// underlying elements
	U *solid.Solid       // u-element
	P *seepage.LiquidGas // pp-element

	// scratchpad. computed @ each ip
	divus float64     // div(us)
	b     []float64   // auxiliary:  b  = a - g = α1・u - ζs - g
	hl    []float64   // auxiliary:  hl = -ρL・b - grad(pl)
	hg    []float64   // temporary:  hg = -ρG・b - grad(pg)
	Kul   [][]float64 // [nu][np] dRus/dpl
	Kug   [][]float64 // [nu][np] dRus/dpl
	Klu   [][]float64 // [np][nu] dRpl/dus
	Kgu   [][]float64 // [np][nu] dRpg/dus

	// for boundary condition derivatives
	Drholdus_ex [][]float64 // [nverts][nverts*ndim] ∂ρl/∂us extrapolted to nodes => if has qlb (flux)
	Drhogdus_ex [][]float64 // [nverts][nverts*ndim] ∂ρg/∂us extrapolted to nodes => if has qgb (flux)
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	ele.SetInfoFunc("solid-liquid-gas", func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *ele.Info {

		// new info
		var info ele.Info

		// u-element info
		u_info := ele.GetInfoFunc("solid")(sim, cell, edat)

		// pp-element info
		edat.Lbb = !sim.Data.NoLBB
		p_info := ele.GetInfoFunc("liquid-gas")(sim, cell, edat)

		// solution variables
		nverts := cell.Shp.Nverts
		info.Dofs = make([][]string, nverts)
		for i, dofs := range u_info.Dofs {
			info.Dofs[i] = append(info.Dofs[i], dofs...)
		}
		for i, dofs := range p_info.Dofs {
			info.Dofs[i] = append(info.Dofs[i], dofs...)
		}

		// maps
		info.Y2F = u_info.Y2F
		for key, val := range p_info.Y2F {
			info.Y2F[key] = val
		}

		// t1 and t2 variables
		info.T1vars = p_info.T1vars
		info.T2vars = u_info.T2vars
		return &info
	})

	// element allocator
	ele.SetAllocator("solid-liquid-gas", func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) ele.Element {

		// basic data
		var o SolidLiquidGas
		o.Sim = sim
		o.Cell = cell
		o.LbbCell = o.Cell
		o.Edat = edat
		o.Ndim = sim.Ndim

		// new LBB cell
		if !sim.Data.NoLBB {
			o.LbbCell = o.Cell.GetSimilar(true)
		}

		// allocate u element
		u_elem := ele.GetAllocator("solid")(sim, cell, edat, x)
		if u_elem == nil {
			chk.Panic("cannot allocate underlying u-element")
		}
		o.U = u_elem.(*solid.Solid)

		// make sure p-element uses the same number of integration points than u-element
		edat.Nip = len(o.U.IpsElem)

		// allocate p-element
		p_elem := ele.GetAllocator("liquid-gas")(sim, o.LbbCell, edat, x)
		if p_elem == nil {
			chk.Panic("cannot allocate underlying p-element")
		}
		o.P = p_elem.(*seepage.LiquidGas)

		// scratchpad. computed @ each ip
		o.b = make([]float64, o.Ndim)
		o.hl = make([]float64, o.Ndim)
		o.hg = make([]float64, o.Ndim)
		o.Kul = la.MatAlloc(o.U.Nu, o.P.Np)
		o.Kug = la.MatAlloc(o.U.Nu, o.P.Np)
		o.Klu = la.MatAlloc(o.P.Np, o.U.Nu)
		o.Kgu = la.MatAlloc(o.P.Np, o.U.Nu)

		// seepage terms
		if o.P.DoExtrap {
			p_nverts := o.P.Cell.Shp.Nverts
			u_nverts := o.U.Cell.Shp.Nverts
			o.Drholdus_ex = la.MatAlloc(p_nverts, u_nverts*o.Ndim)
			o.Drhogdus_ex = la.MatAlloc(p_nverts, u_nverts*o.Ndim)
		}

		// return new element
		return &o
	})
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o *SolidLiquidGas) Id() int { return o.Cell.Id }

// SetEqs set equations
func (o *SolidLiquidGas) SetEqs(eqs [][]int, mixedform_eqs []int) (err error) {

	// u: equations
	u_info := ele.GetInfoFunc("solid")(o.Sim, o.Cell, o.Edat)
	u_nverts := len(u_info.Dofs)
	u_eqs := make([][]int, u_nverts)
	for i := 0; i < u_nverts; i++ {
		nkeys := len(u_info.Dofs[i])
		u_eqs[i] = make([]int, nkeys)
		for j := 0; j < nkeys; j++ {
			u_eqs[i][j] = eqs[i][j]
		}
	}

	// p: equations
	p_info := ele.GetInfoFunc("liquid-gas")(o.Sim, o.LbbCell, o.Edat)
	p_nverts := len(p_info.Dofs)
	p_eqs := make([][]int, p_nverts)
	for i := 0; i < p_nverts; i++ {
		start := len(u_info.Dofs[i])
		nkeys := len(p_info.Dofs[i])
		p_eqs[i] = make([]int, nkeys)
		for j := 0; j < nkeys; j++ {
			p_eqs[i][j] = eqs[i][start+j]
		}
	}

	// set equations
	err = o.U.SetEqs(u_eqs, mixedform_eqs)
	if err != nil {
		return
	}
	return o.P.SetEqs(p_eqs, nil)
}

// SetEleConds set element conditions
func (o *SolidLiquidGas) SetEleConds(key string, f dbf.T, extra string) (err error) {
	err = o.U.SetEleConds(key, f, extra)
	if err != nil {
		return
	}
	return o.P.SetEleConds(key, f, extra)
}

// InterpStarVars interpolates star variables to integration points
func (o *SolidLiquidGas) InterpStarVars(sol *ele.Solution) (err error) {
	err = o.U.InterpStarVars(sol)
	if err != nil {
		return
	}
	return o.P.InterpStarVars(sol)
}

// adds -R to global residual vector fb
func (o *SolidLiquidGas) AddToRhs(fb []float64, sol *ele.Solution) (err error) {

	// clear variables
	if o.P.DoExtrap {
		la.VecFill(o.P.Rhol_ex, 0)
		la.VecFill(o.P.Rhog_ex, 0)
	}
	if o.U.UseB {
		la.VecFill(o.U.Fi, 0)
	}

	// for each integration point
	O := o.P.LgsVars
	α4 := sol.DynCfs.GetAlp4()
	β1 := sol.DynCfs.GetBet1()
	u_nverts := o.U.Cell.Shp.Nverts
	p_nverts := o.P.Cell.Shp.Nverts
	var coef, plt, pgt, klr, kgr, ρl, ρg, ρ, p, divvs float64
	for idx, ip := range o.U.IpsElem {

		// interpolation functions, gradients and variables @ ip
		err = o.ipvars(idx, sol)
		if err != nil {
			return
		}
		coef = o.U.Cell.Shp.J * ip[3]
		S := o.U.Cell.Shp.S
		G := o.U.Cell.Shp.G
		Sb := o.P.Cell.Shp.S
		Gb := o.P.Cell.Shp.G

		// axisymmetric case
		radius := 1.0
		if sol.Axisym {
			radius = o.U.Cell.Shp.AxisymGetRadius(o.U.X)
			coef *= radius
		}

		// tpm variables
		plt = β1*o.P.Pl - o.P.PsiL[idx]
		pgt = β1*o.P.Pg - o.P.PsiG[idx]
		klr = o.P.Mdl.Cnd.Klr(o.P.States[idx].A_sl)
		kgr = o.P.Mdl.Cnd.Kgr(1.0 - o.P.States[idx].A_sl)
		err = o.P.Mdl.CalcLgs(O, o.P.States[idx], o.P.Pl, o.P.Pg, o.divus, false)
		if err != nil {
			return
		}
		ρl = O.A_ρl
		ρg = O.A_ρg
		ρ = O.A_ρ
		p = O.A_p
		σe := o.U.States[idx].Sig
		divvs = α4*o.divus - o.U.DivChi[idx]

		// compute augmented filter velocities
		for i := 0; i < o.Ndim; i++ {
			o.P.Wlb[i], o.P.Wgb[i] = 0, 0
			for j := 0; j < o.Ndim; j++ {
				o.P.Wlb[i] += klr * o.P.Mdl.Klsat[i][j] * o.hl[j]
				o.P.Wgb[i] += kgr * o.P.Mdl.Kgsat[i][j] * o.hg[j]
			}
		}

		// p: add negative of residual term to fb
		for m := 0; m < p_nverts; m++ {
			rl := o.P.Plmap[m]
			rg := o.P.Pgmap[m]
			fb[rl] -= coef * Sb[m] * (O.Cpl*plt + O.Cpg*pgt + O.Cvs*divvs)
			fb[rg] -= coef * Sb[m] * (O.Dpl*plt + O.Dpg*pgt + O.Dvs*divvs)
			for i := 0; i < o.Ndim; i++ {
				fb[rl] += coef * Gb[m][i] * o.P.Wlb[i] // += coef * div(ρl*wl)
				fb[rg] += coef * Gb[m][i] * o.P.Wgb[i] // += coef * div(ρg*wg)
			}
			if o.P.DoExtrap {
				o.P.Rhol_ex[m] += o.P.Emat[m][idx] * ρl
				o.P.Rhog_ex[m] += o.P.Emat[m][idx] * ρg
			}
		}

		// u: add negative of residual term to fb; see Eqs. (38b) and (45b) [1]
		if o.U.UseB {
			solid.IpBmatrix(o.U.B, o.Ndim, u_nverts, G, radius, S, sol.Axisym)
			la.MatTrVecMulAdd(o.U.Fi, coef, o.U.B, σe) // fi += coef * tr(B) * σ
			for m := 0; m < u_nverts; m++ {
				for i := 0; i < o.Ndim; i++ {
					ru := o.U.Umap[i+m*o.Ndim]
					fb[ru] -= coef * S[m] * ρ * o.b[i]
					fb[ru] += coef * p * G[m][i]
				}
			}
		} else {
			for m := 0; m < u_nverts; m++ {
				for i := 0; i < o.Ndim; i++ {
					ru := o.U.Umap[i+m*o.Ndim]
					fb[ru] -= coef * S[m] * ρ * o.b[i]
					for j := 0; j < o.Ndim; j++ {
						fb[ru] -= coef * tsr.M2T(σe, i, j) * G[m][j]
					}
					fb[ru] += coef * p * G[m][i]
				}
			}
		}
	}

	// add fi term to fb, if using B matrix
	if o.U.UseB {
		for i, I := range o.U.Umap {
			fb[I] -= o.U.Fi[i]
		}
	}

	// external forces
	if len(o.U.NatBcs) > 0 {
		err = o.U.AddSurfLoadsToRhs(fb, sol)
		if err != nil {
			return
		}
	}

	// contribution from natural boundary conditions
	if len(o.P.NatBcs) > 0 {
		return o.P.AddNatBcsToRhs(fb, sol)
	}
	return
}

// adds element K to global Jacobian matrix Kb
func (o *SolidLiquidGas) AddToKb(Kb *la.Triplet, sol *ele.Solution, firstIt bool) (err error) {

	// clear matrices
	u_nverts := o.U.Cell.Shp.Nverts
	p_nverts := o.P.Cell.Shp.Nverts
	for i := 0; i < o.P.Np; i++ {
		for j := 0; j < o.P.Np; j++ {
			o.P.Kll[i][j], o.P.Klg[i][j] = 0, 0
			o.P.Kgl[i][j], o.P.Kgg[i][j] = 0, 0
		}
	}
	for i := 0; i < o.U.Nu; i++ {
		for j := 0; j < o.P.Np; j++ {
			o.Kul[i][j], o.Kug[i][j] = 0, 0
			o.Klu[j][i], o.Kgu[j][i] = 0, 0
		}
		for j := 0; j < o.U.Nu; j++ {
			o.U.K[i][j] = 0
		}
	}
	if o.P.DoExtrap {
		for i := 0; i < p_nverts; i++ {
			o.P.Rhol_ex[i] = 0
			for j := 0; j < p_nverts; j++ {
				o.P.Drholdpl_ex[i][j] = 0
				o.P.Drhogdpg_ex[i][j] = 0
			}
			for j := 0; j < o.U.Nu; j++ {
				o.Drholdus_ex[i][j] = 0
				o.Drhogdus_ex[i][j] = 0
			}
		}
	}

	// for each integration point
	O := o.P.LgsVars
	Cl := o.P.Mdl.Liq.C
	Cg := o.P.Mdl.Gas.C
	α1 := sol.DynCfs.GetAlp1()
	α4 := sol.DynCfs.GetAlp4()
	β1 := sol.DynCfs.GetBet1()
	var coef, plt, pgt, klr, kgr, ρL, ρG, ρl, ρg, ρ, divvs, dhldpl_nj, dhgdpg_nj float64
	for idx, ip := range o.U.IpsElem {

		// interpolation functions, gradients and variables @ ip
		err = o.ipvars(idx, sol)
		if err != nil {
			return
		}
		coef = o.U.Cell.Shp.J * ip[3]
		S := o.U.Cell.Shp.S
		G := o.U.Cell.Shp.G
		Sb := o.P.Cell.Shp.S
		Gb := o.P.Cell.Shp.G

		// axisymmetric case
		radius := 1.0
		if sol.Axisym {
			radius = o.U.Cell.Shp.AxisymGetRadius(o.U.X)
			coef *= radius
		}

		// tpm variables
		plt = β1*o.P.Pl - o.P.PsiL[idx]
		pgt = β1*o.P.Pg - o.P.PsiG[idx]
		klr = o.P.Mdl.Cnd.Klr(o.P.States[idx].A_sl)
		kgr = o.P.Mdl.Cnd.Kgr(1.0 - o.P.States[idx].A_sl)
		ρL = o.P.States[idx].A_ρL
		ρG = o.P.States[idx].A_ρG
		err = o.P.Mdl.CalcLgs(O, o.P.States[idx], o.P.Pl, o.P.Pg, o.divus, true)
		if err != nil {
			return
		}
		ρl = O.A_ρl
		ρg = O.A_ρg
		ρ = O.A_ρ
		divvs = α4*o.divus - o.U.DivChi[idx]

		// Jacobian
		for n := 0; n < p_nverts; n++ {

			// clear auxiliary variables
			for i := 0; i < o.Ndim; i++ {
				o.P.Dwlbdpl_n[i], o.P.Dwlbdpg_n[i] = 0, 0
				o.P.Dwgbdpl_n[i], o.P.Dwgbdpg_n[i] = 0, 0
			}

			// loop over space dimension
			for j := 0; j < o.Ndim; j++ {

				// {pl,pg,us} versus {pl,pg,us}
				for m := 0; m < u_nverts; m++ {
					c := j + m*o.Ndim

					// ∂rβ/∂u^m
					o.Klu[n][c] += coef * Sb[n] * (O.DCplduM*plt + O.DCpgduM*pgt + α4*O.Cvs) * G[m][j]
					o.Kgu[n][c] += coef * Sb[n] * (O.DDplduM*plt + O.DDpgduM*pgt + α4*O.Dvs) * G[m][j]

					// ∂wβb/∂u^m
					for i := 0; i < o.Ndim; i++ {
						o.Klu[n][c] += coef * Gb[n][i] * S[m] * α1 * ρL * klr * o.P.Mdl.Klsat[i][j]
						o.Kgu[n][c] += coef * Gb[n][i] * S[m] * α1 * ρG * kgr * o.P.Mdl.Kgsat[i][j]
					}

					// ∂ru/∂pβ^n  and  ∂p/∂pβ^n
					o.Kul[c][n] += coef * (S[m]*Sb[n]*O.Dρdpl*o.b[j] - G[m][j]*Sb[n]*O.Dpdpl)
					o.Kug[c][n] += coef * (S[m]*Sb[n]*O.Dρdpg*o.b[j] - G[m][j]*Sb[n]*O.Dpdpg)

					// extrapolation term
					if o.P.DoExtrap {
						o.Drholdus_ex[n][c] += o.P.Emat[n][idx] * O.DρlduM * G[m][j]
						o.Drhogdus_ex[n][c] += o.P.Emat[n][idx] * O.DρgduM * G[m][j]
					}
				}

				// compute auxiliary derivatives
				dhldpl_nj = Sb[n]*Cl*o.P.Grav[j] - Gb[n][j]
				dhgdpg_nj = Sb[n]*Cg*o.P.Grav[j] - Gb[n][j]
				for i := 0; i < o.Ndim; i++ {
					o.P.Dwlbdpl_n[i] += o.P.Mdl.Klsat[i][j] * (Sb[n]*O.Dklrdpl*o.hl[j] + klr*dhldpl_nj)
					o.P.Dwlbdpg_n[i] += o.P.Mdl.Klsat[i][j] * (Sb[n] * O.Dklrdpg * o.hl[j])
					o.P.Dwgbdpl_n[i] += o.P.Mdl.Kgsat[i][j] * (Sb[n] * O.Dkgrdpl * o.hg[j])
					o.P.Dwgbdpg_n[i] += o.P.Mdl.Kgsat[i][j] * (Sb[n]*O.Dkgrdpg*o.hg[j] + kgr*dhgdpg_nj)
				}
			}

			// {pl,pg} versus {pl,pg}
			for m := 0; m < p_nverts; m++ {

				// ∂rβ/∂pγ
				o.P.Kll[m][n] += coef * Sb[m] * Sb[n] * (O.DCpldpl*plt + O.DCpgdpl*pgt + O.DCvsdpl*divvs + β1*O.Cpl)
				o.P.Klg[m][n] += coef * Sb[m] * Sb[n] * (O.DCpldpg*plt + O.DCpgdpg*pgt + O.DCvsdpg*divvs + β1*O.Cpg)
				o.P.Kgl[m][n] += coef * Sb[m] * Sb[n] * (O.DDpldpl*plt + O.DDpgdpl*pgt + O.DDvsdpl*divvs + β1*O.Dpl)
				o.P.Kgg[m][n] += coef * Sb[m] * Sb[n] * (O.DDpldpg*plt + O.DDpgdpg*pgt + O.DDvsdpg*divvs + β1*O.Dpg)

				// ∂wβb/∂pγ
				for i := 0; i < o.Ndim; i++ {
					o.P.Kll[m][n] -= coef * Gb[m][i] * o.P.Dwlbdpl_n[i]
					o.P.Klg[m][n] -= coef * Gb[m][i] * o.P.Dwlbdpg_n[i]
					o.P.Kgl[m][n] -= coef * Gb[m][i] * o.P.Dwgbdpl_n[i]
					o.P.Kgg[m][n] -= coef * Gb[m][i] * o.P.Dwgbdpg_n[i]
				}

				// extrapolation term
				if o.P.DoExtrap {
					o.P.Drholdpl_ex[m][n] += o.P.Emat[m][idx] * O.Cpl * Sb[n]
					o.P.Drhogdpg_ex[m][n] += o.P.Emat[m][idx] * O.Dpg * Sb[n]
				}
			}

			// extrapolation terms
			if o.P.DoExtrap {
				o.P.Rhol_ex[n] += o.P.Emat[n][idx] * ρl
				o.P.Rhog_ex[n] += o.P.Emat[n][idx] * ρg
			}
		}

		// {u} versus {u}
		for m := 0; m < u_nverts; m++ {
			for i := 0; i < o.Ndim; i++ {
				r := i + m*o.Ndim
				for n := 0; n < u_nverts; n++ {
					for j := 0; j < o.Ndim; j++ {
						c := j + n*o.Ndim
						o.U.K[r][c] += coef * S[m] * (S[n]*α1*ρ*tsr.It[i][j] + O.DρduM*o.b[i]*G[n][j])
					}
				}
			}
		}

		// consistent tangent model matrix
		err = o.U.MdlSmall.CalcD(o.U.D, o.U.States[idx], firstIt)
		if err != nil {
			return
		}

		// Kuu: add stiffness term ∂(σe・G^m)/∂us^n
		if o.U.UseB {
			solid.IpBmatrix(o.U.B, o.Ndim, u_nverts, G, radius, S, sol.Axisym)
			la.MatTrMulAdd3(o.U.K, coef, o.U.B, o.U.D, o.U.B) // K += coef * tr(B) * D * B
		} else {
			solid.IpAddToKt(o.U.K, u_nverts, o.Ndim, coef, G, o.U.D)
		}
	}

	// contribution from natural boundary conditions
	if o.P.HasSeep {
		err = o.P.AddNatBcsToJac(sol)
		if err != nil {
			return
		}
		err = o.AddNatBcsToJac(sol)
		if err != nil {
			return
		}
	}

	// assemble K matrices into Kb
	o.P.AssembleKs(Kb)
	for i, I := range o.P.Plmap {
		for j, J := range o.U.Umap {
			Kb.Put(I, J, o.Klu[i][j])
			Kb.Put(J, I, o.Kul[j][i])
		}
	}
	for i, I := range o.P.Pgmap {
		for j, J := range o.U.Umap {
			Kb.Put(I, J, o.Kgu[i][j])
			Kb.Put(J, I, o.Kug[j][i])
		}
	}
	for i, I := range o.U.Umap {
		for j, J := range o.U.Umap {
			Kb.Put(I, J, o.U.K[i][j])
		}
	}
	return
}

// Update perform (tangent) update
func (o *SolidLiquidGas) Update(sol *ele.Solution) (err error) {
	err = o.U.Update(sol)
	if err != nil {
		return
	}
	return o.P.Update(sol)
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// SetIniIvs sets initial ivs for given values in sol and ivs map
func (o *SolidLiquidGas) SetIniIvs(sol *ele.Solution, ivs map[string][]float64) (err error) {

	// set p-element first
	err = o.P.SetIniIvs(sol, nil)
	if err != nil {
		return
	}

	// initial stresses given
	if _, okk := ivs["svT"]; okk {

		// total vertical stresses and K0
		nip := len(o.U.IpsElem)
		svT := ivs["svT"]
		K0s := ivs["K0"]
		chk.IntAssert(len(svT), nip)
		chk.IntAssert(len(K0s), 1)
		K0 := K0s[0]

		// for each integration point
		sx := make([]float64, nip)
		sy := make([]float64, nip)
		sz := make([]float64, nip)
		for i, ip := range o.U.IpsElem {

			// compute pl and pg @ ip
			err = o.P.Cell.Shp.CalcAtIp(o.P.X, ip, false)
			if err != nil {
				return
			}
			pl, pg := 0.0, 0.0
			for m := 0; m < o.P.Cell.Shp.Nverts; m++ {
				rl := o.P.Plmap[m]
				rg := o.P.Pgmap[m]
				pl += o.P.Cell.Shp.S[m] * sol.Y[rl]
				pg += o.P.Cell.Shp.S[m] * sol.Y[rg]
			}

			// compute effective stresses
			sl := o.P.States[i].A_sl
			sg := 1.0 - sl
			p := pl*sl + pg*sg
			svE := svT[i] + p
			shE := K0 * svE
			sx[i], sy[i], sz[i] = shE, svE, shE
			if o.Ndim == 3 {
				sx[i], sy[i], sz[i] = shE, shE, svE
			}
		}
		ivs = map[string][]float64{"sx": sx, "sy": sy, "sz": sz}
	}

	// set u-element
	return o.U.SetIniIvs(sol, ivs)
}

// BackupIvs create copy of internal variables
func (o *SolidLiquidGas) BackupIvs(aux bool) (err error) {
	err = o.U.BackupIvs(aux)
	if err != nil {
		return
	}
	return o.P.BackupIvs(aux)
}

// RestoreIvs restore internal variables from copies
func (o *SolidLiquidGas) RestoreIvs(aux bool) (err error) {
	err = o.U.RestoreIvs(aux)
	if err != nil {
		return
	}
	return o.P.RestoreIvs(aux)
}

// Ureset fixes internal variables after u (displacements) have been zeroed
func (o *SolidLiquidGas) Ureset(sol *ele.Solution) (err error) {
	u_nverts := o.U.Cell.Shp.Nverts
	for idx, ip := range o.U.IpsElem {
		err = o.U.Cell.Shp.CalcAtIp(o.U.X, ip, true)
		if err != nil {
			return
		}
		G := o.U.Cell.Shp.G
		var divus float64
		for m := 0; m < u_nverts; m++ {
			for i := 0; i < o.Ndim; i++ {
				r := o.U.Umap[i+m*o.Ndim]
				divus += G[m][i] * sol.Y[r]
			}
		}
		o.P.States[idx].A_ns0 = (1.0 - divus) * (1.0 - o.P.Mdl.Nf0)
		o.P.StatesBkp[idx].A_ns0 = o.P.States[idx].A_ns0
	}
	err = o.U.Ureset(sol)
	if err != nil {
		return
	}
	return o.P.Ureset(sol)
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o *SolidLiquidGas) Encode(enc utl.Encoder) (err error) {
	err = o.U.Encode(enc)
	if err != nil {
		return
	}
	return o.P.Encode(enc)
}

// Decode decodes internal variables
func (o *SolidLiquidGas) Decode(dec utl.Decoder) (err error) {
	err = o.U.Decode(dec)
	if err != nil {
		return
	}
	return o.P.Decode(dec)
}

// OutIpCoords returns the coordinates of integration points
func (o *SolidLiquidGas) OutIpCoords() (C [][]float64) {
	return o.U.OutIpCoords()
}

// OutIpKeys returns the integration points' keys
func (o *SolidLiquidGas) OutIpKeys() []string {
	keys := append(o.U.OutIpKeys(), "nf", "pl", "pg", "pc", "sl", "RhoL", "RhoG")
	keys = append(keys, seepage.LiqFlowKeys(o.Ndim)...)
	return append(keys, seepage.GasFlowKeys(o.Ndim)...)
}

// OutIpVals returns the integration points' values corresponding to keys
func (o *SolidLiquidGas) OutIpVals(M *ele.IpsMap, sol *ele.Solution) {
	o.U.OutIpVals(M, sol)
	flowL := seepage.LiqFlowKeys(o.Ndim)
	flowG := seepage.GasFlowKeys(o.Ndim)
	nip := len(o.U.IpsElem)
	for idx, _ := range o.U.IpsElem {
		err := o.ipvars(idx, sol)
		if err != nil {
			return
		}
		ns := (1.0 - o.divus) * o.P.States[idx].A_ns0
		sl := o.P.States[idx].A_sl
		sg := 1.0 - sl
		ρL := o.P.States[idx].A_ρL
		ρG := o.P.States[idx].A_ρG
		klr := o.P.Mdl.Cnd.Klr(sl)
		kgr := o.P.Mdl.Cnd.Klr(sg)
		M.Set("nf", idx, nip, 1.0-ns)
		M.Set("pl", idx, nip, o.P.Pl)
		M.Set("pg", idx, nip, o.P.Pg)
		M.Set("pc", idx, nip, o.P.Pg-o.P.Pl)
		M.Set("sl", idx, nip, sl)
		M.Set("RhoL", idx, nip, ρL)
		M.Set("RhoG", idx, nip, ρG)
		for i := 0; i < o.Ndim; i++ {
			var nwl_i, nwg_i float64
			for j := 0; j < o.Ndim; j++ {
				nwl_i += klr * o.P.Mdl.Klsat[i][j] * (o.P.Grav[j] - o.P.GradPl[j]/ρL)
				nwg_i += kgr * o.P.Mdl.Kgsat[i][j] * (o.P.Grav[j] - o.P.GradPg[j]/ρG)
			}
			M.Set(flowL[i], idx, nip, nwl_i)
			M.Set(flowG[i], idx, nip, nwg_i)
		}
	}
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// ipvars computes current values @ integration points. idx == index of integration point
func (o *SolidLiquidGas) ipvars(idx int, sol *ele.Solution) (err error) {

	// interpolation functions and gradients
	err = o.P.Cell.Shp.CalcAtIp(o.P.X, o.U.IpsElem[idx], true)
	if err != nil {
		return
	}
	err = o.U.Cell.Shp.CalcAtIp(o.U.X, o.U.IpsElem[idx], true)
	if err != nil {
		return
	}

	// auxiliary
	ρL := o.P.States[idx].A_ρL
	ρG := o.P.States[idx].A_ρG
	o.P.ComputeGrav(sol.T)

	// clear gpl and recover u-variables @ ip
	o.divus = 0
	for i := 0; i < o.Ndim; i++ {
		o.P.GradPl[i], o.P.GradPg[i], o.U.Us[i] = 0, 0, 0
		for m := 0; m < o.U.Cell.Shp.Nverts; m++ {
			r := o.U.Umap[i+m*o.Ndim]
			o.U.Us[i] += o.U.Cell.Shp.S[m] * sol.Y[r]
			o.divus += o.U.Cell.Shp.G[m][i] * sol.Y[r]
		}
	}

	// recover p-variables @ ip
	o.P.Pl, o.P.Pg = 0, 0
	for m := 0; m < o.P.Cell.Shp.Nverts; m++ {
		rl := o.P.Plmap[m]
		rg := o.P.Pgmap[m]
		o.P.Pl += o.P.Cell.Shp.S[m] * sol.Y[rl]
		o.P.Pg += o.P.Cell.Shp.S[m] * sol.Y[rg]
		for i := 0; i < o.Ndim; i++ {
			o.P.GradPl[i] += o.P.Cell.Shp.G[m][i] * sol.Y[rl]
			o.P.GradPg[i] += o.P.Cell.Shp.G[m][i] * sol.Y[rg]
		}
	}

	// compute b, hl and hg
	α1 := sol.DynCfs.GetAlp1()
	for i := 0; i < o.Ndim; i++ {
		o.b[i] = α1*o.U.Us[i] - o.U.Zet[idx][i] - o.P.Grav[i]
		o.hl[i] = -ρL*o.b[i] - o.P.GradPl[i]
		o.hg[i] = -ρG*o.b[i] - o.P.GradPg[i]
	}
	return
}

// AddNatBcsToJac adds contribution from natural boundary conditions to Jacobian
func (o *SolidLiquidGas) AddNatBcsToJac(sol *ele.Solution) (err error) {

	// compute surface integral
	u_nverts := o.U.Cell.Shp.Nverts
	var shift float64
	var pl, fl, plmax, g, rmp float64
	for idx, nbc := range o.P.NatBcs {

		// plmax shift
		shift = nbc.Fcn.F(sol.T, nil)

		// loop over ips of face
		for jdx, ipf := range o.P.IpsFace {

			// interpolation functions and gradients @ face
			iface := nbc.IdxFace
			err = o.P.Cell.Shp.CalcAtFaceIp(o.P.X, ipf, iface)
			if err != nil {
				return
			}
			Sf := o.P.Cell.Shp.Sf
			Jf := la.VecNorm(o.P.Cell.Shp.Fnvec)
			coef := ipf[3] * Jf

			// select natural boundary condition type
			switch nbc.Key {
			case "seep":

				// variables extrapolated to face
				_, pl, fl = o.P.CalcFaceIpVars(iface, sol)
				plmax = o.P.Plmax[idx][jdx] - shift
				if plmax < 0 {
					plmax = 0
				}

				// compute derivatives
				g = pl - plmax // Eq. (24)
				rmp = o.P.Ramp(fl + o.P.Kap*g)
				for i, m := range o.P.Cell.Shp.FaceLocalVerts[iface] {
					for n := 0; n < u_nverts; n++ {
						for j := 0; j < o.Ndim; j++ {
							c := j + n*o.Ndim
							for l, r := range o.P.Cell.Shp.FaceLocalVerts[iface] {
								o.Klu[m][c] += coef * Sf[i] * Sf[l] * o.Drholdus_ex[r][c] * rmp
							}
						}
					}
				}
			}
		}
	}
	return
}
