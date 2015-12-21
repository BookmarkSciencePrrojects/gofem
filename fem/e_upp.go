// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/tsr"
)

// ElemUPP implements the u-pl-pg formulation
type ElemUPP struct {

	// auxiliary
	Sim     *inp.Simulation // simulation
	Cell    *inp.Cell       // cell
	LbbCell *inp.Cell       // if LBB==false, same as Cell; otherwise LbbCell is a new cell with less vertices
	Edat    *inp.ElemData   // element data; stored in allocator to be used in Connect
	Ndim    int             // space dimension

	// underlying elements
	U *ElemU  // u-element
	P *ElemPP // pp-element

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
	dρldus_ex [][]float64 // [nverts][nverts*ndim] ∂ρl/∂us extrapolted to nodes => if has qlb (flux)
	dρgdus_ex [][]float64 // [nverts][nverts*ndim] ∂ρg/∂us extrapolted to nodes => if has qgb (flux)
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	infogetters["upp"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *Info {

		// new info
		var info Info

		// u-element info
		u_info := infogetters["u"](sim, cell, edat)

		// p-element info
		edat.Lbb = !sim.Data.NoLBB
		p_info := infogetters["pp"](sim, cell, edat)

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
	}

	// element allocator
	eallocators["upp"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) Elem {

		// basic data
		var o ElemUPP
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
		u_elem := eallocators["u"](sim, cell, edat, x)
		if u_elem == nil {
			chk.Panic("cannot allocate underlying u-element")
		}
		o.U = u_elem.(*ElemU)

		// make sure p-element uses the same number of integration points than u-element
		edat.Nip = len(o.U.IpsElem)

		// allocate p-element
		p_elem := eallocators["pp"](sim, o.LbbCell, edat, x)
		if p_elem == nil {
			chk.Panic("cannot allocate underlying p-element")
		}
		o.P = p_elem.(*ElemPP)

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
			o.dρldus_ex = la.MatAlloc(p_nverts, u_nverts*o.Ndim)
			o.dρgdus_ex = la.MatAlloc(p_nverts, u_nverts*o.Ndim)
		}

		// return new element
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o *ElemUPP) Id() int { return o.Cell.Id }

// SetEqs set equations
func (o *ElemUPP) SetEqs(eqs [][]int, mixedform_eqs []int) (err error) {

	// u: equations
	u_info := infogetters["u"](o.Sim, o.Cell, o.Edat)
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
	p_info := infogetters["pp"](o.Sim, o.LbbCell, o.Edat)
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
func (o *ElemUPP) SetEleConds(key string, f fun.Func, extra string) (err error) {
	err = o.U.SetEleConds(key, f, extra)
	if err != nil {
		return
	}
	return o.P.SetEleConds(key, f, extra)
}

// InterpStarVars interpolates star variables to integration points
func (o *ElemUPP) InterpStarVars(sol *Solution) (err error) {

	// for each integration point
	u_nverts := o.U.Cell.Shp.Nverts
	p_nverts := o.P.Cell.Shp.Nverts
	for idx, ip := range o.U.IpsElem {

		// interpolation functions and gradients
		err = o.P.Cell.Shp.CalcAtIp(o.P.X, ip, true)
		if err != nil {
			return
		}
		err = o.U.Cell.Shp.CalcAtIp(o.U.X, ip, true)
		if err != nil {
			return
		}
		S := o.U.Cell.Shp.S
		G := o.U.Cell.Shp.G
		Sb := o.P.Cell.Shp.S

		// clear local variables
		o.P.ψl[idx], o.P.ψg[idx], o.U.divχs[idx] = 0, 0, 0
		for i := 0; i < o.Ndim; i++ {
			o.U.ζs[idx][i], o.U.χs[idx][i] = 0, 0
		}

		// p-variables
		for m := 0; m < p_nverts; m++ {
			rl := o.P.Plmap[m]
			rg := o.P.Pgmap[m]
			o.P.ψl[idx] += Sb[m] * sol.Psi[rl]
			o.P.ψg[idx] += Sb[m] * sol.Psi[rg]
		}

		// u-variables
		for m := 0; m < u_nverts; m++ {
			for i := 0; i < o.Ndim; i++ {
				r := o.U.Umap[i+m*o.Ndim]
				o.U.ζs[idx][i] += S[m] * sol.Zet[r]
				o.U.χs[idx][i] += S[m] * sol.Chi[r]
				o.U.divχs[idx] += G[m][i] * sol.Chi[r]
			}
		}
	}
	return
}

// adds -R to global residual vector fb
func (o *ElemUPP) AddToRhs(fb []float64, sol *Solution) (err error) {

	// clear variables
	if o.P.DoExtrap {
		la.VecFill(o.P.ρl_ex, 0)
		la.VecFill(o.P.ρg_ex, 0)
	}
	if o.U.UseB {
		la.VecFill(o.U.fi, 0)
	}

	// for each integration point
	O := o.P.res
	α4 := sol.DynCfs.α4
	β1 := sol.DynCfs.β1
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
		plt = β1*o.P.pl - o.P.ψl[idx]
		pgt = β1*o.P.pg - o.P.ψg[idx]
		klr = o.P.Mdl.Cnd.Klr(o.P.States[idx].A_sl)
		kgr = o.P.Mdl.Cnd.Kgr(1.0 - o.P.States[idx].A_sl)
		err = o.P.Mdl.CalcLgs(o.P.res, o.P.States[idx], o.P.pl, o.P.pg, o.divus, false)
		if err != nil {
			return
		}
		ρl = o.P.res.A_ρl
		ρg = o.P.res.A_ρg
		ρ = o.P.res.A_ρ
		p = o.P.res.A_p
		σe := o.U.States[idx].Sig
		divvs = α4*o.divus - o.U.divχs[idx]

		// compute augmented filter velocities
		for i := 0; i < o.Ndim; i++ {
			o.P.wlb[i], o.P.wgb[i] = 0, 0
			for j := 0; j < o.Ndim; j++ {
				o.P.wlb[i] += klr * o.P.Mdl.Klsat[i][j] * o.hl[j]
				o.P.wgb[i] += kgr * o.P.Mdl.Kgsat[i][j] * o.hg[j]
			}
		}

		// p: add negative of residual term to fb
		for m := 0; m < p_nverts; m++ {
			rl := o.P.Plmap[m]
			rg := o.P.Pgmap[m]
			fb[rl] -= coef * Sb[m] * (O.Cpl*plt + O.Cpg*pgt + O.Cvs*divvs)
			fb[rg] -= coef * Sb[m] * (O.Dpl*plt + O.Dpg*pgt + O.Dvs*divvs)
			for i := 0; i < o.Ndim; i++ {
				fb[rl] += coef * Gb[m][i] * o.P.wlb[i] // += coef * div(ρl*wl)
				fb[rg] += coef * Gb[m][i] * o.P.wgb[i] // += coef * div(ρg*wg)
			}
			if o.P.DoExtrap {
				o.P.ρl_ex[m] += o.P.Emat[m][idx] * ρl
				o.P.ρg_ex[m] += o.P.Emat[m][idx] * ρg
			}
		}

		// u: add negative of residual term to fb; see Eqs. (38b) and (45b) [1]
		if o.U.UseB {
			IpBmatrix(o.U.B, o.Ndim, u_nverts, G, radius, S, sol.Axisym)
			la.MatTrVecMulAdd(o.U.fi, coef, o.U.B, σe) // fi += coef * tr(B) * σ
			for m := 0; m < u_nverts; m++ {
				for i := 0; i < o.Ndim; i++ {
					r := o.U.Umap[i+m*o.Ndim]
					fb[r] -= coef * S[m] * ρ * o.b[i]
					fb[r] += coef * p * G[m][i]
				}
			}
		} else {
			for m := 0; m < u_nverts; m++ {
				for i := 0; i < o.Ndim; i++ {
					r := o.U.Umap[i+m*o.Ndim]
					fb[r] -= coef * S[m] * ρ * o.b[i]
					for j := 0; j < o.Ndim; j++ {
						fb[r] -= coef * tsr.M2T(σe, i, j) * G[m][j]
					}
					fb[r] += coef * p * G[m][i]
				}
			}
		}
	}

	// add fi term to fb, if using B matrix
	if o.U.UseB {
		for i, I := range o.U.Umap {
			fb[I] -= o.U.fi[i]
		}
	}

	// external forces
	if len(o.U.NatBcs) > 0 {
		err = o.U.add_surfloads_to_rhs(fb, sol)
		if err != nil {
			return
		}
	}

	// contribution from natural boundary conditions
	if len(o.P.NatBcs) > 0 {
		return o.P.add_natbcs_to_rhs(fb, sol)
	}
	return
}

// adds element K to global Jacobian matrix Kb
func (o *ElemUPP) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (err error) {

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
			o.P.ρl_ex[i] = 0
			for j := 0; j < p_nverts; j++ {
				o.P.dρldpl_ex[i][j] = 0
				o.P.dρgdpg_ex[i][j] = 0
			}
			for j := 0; j < o.U.Nu; j++ {
				o.dρldus_ex[i][j] = 0
				o.dρgdus_ex[i][j] = 0
			}
		}
	}

	// for each integration point
	O := o.P.res
	Cl := o.P.Mdl.Liq.C
	Cg := o.P.Mdl.Gas.C
	α1 := sol.DynCfs.α1
	α4 := sol.DynCfs.α4
	β1 := sol.DynCfs.β1
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
		plt = β1*o.P.pl - o.P.ψl[idx]
		pgt = β1*o.P.pg - o.P.ψg[idx]
		klr = o.P.Mdl.Cnd.Klr(o.P.States[idx].A_sl)
		kgr = o.P.Mdl.Cnd.Kgr(1.0 - o.P.States[idx].A_sl)
		ρL = o.P.States[idx].A_ρL
		ρG = o.P.States[idx].A_ρG
		err = o.P.Mdl.CalcLgs(o.P.res, o.P.States[idx], o.P.pl, o.P.pg, o.divus, false)
		if err != nil {
			return
		}
		ρl = o.P.res.A_ρl
		ρg = o.P.res.A_ρg
		ρ = o.P.res.A_ρ
		divvs = α4*o.divus - o.U.divχs[idx]

		// Jacobian
		for n := 0; n < p_nverts; n++ {

			// clear auxiliary variables
			for i := 0; i < o.Ndim; i++ {
				o.P.dwlbdpl_n[i], o.P.dwlbdpg_n[i] = 0, 0
				o.P.dwgbdpl_n[i], o.P.dwgbdpg_n[i] = 0, 0
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
						o.dρldus_ex[n][c] += o.P.Emat[n][idx] * O.DρlduM * G[m][j]
						o.dρgdus_ex[n][c] += o.P.Emat[n][idx] * O.DρgduM * G[m][j]
					}
				}

				// compute auxiliary derivatives
				dhldpl_nj = Sb[n]*Cl*o.P.g[j] - Gb[n][j]
				dhgdpg_nj = Sb[n]*Cg*o.P.g[j] - Gb[n][j]
				for i := 0; i < o.Ndim; i++ {
					o.P.dwlbdpl_n[i] += o.P.Mdl.Klsat[i][j] * (Sb[n]*O.Dklrdpl*o.hl[j] + klr*dhldpl_nj)
					o.P.dwlbdpg_n[i] += o.P.Mdl.Klsat[i][j] * (Sb[n] * O.Dklrdpg * o.hl[j])
					o.P.dwgbdpl_n[i] += o.P.Mdl.Kgsat[i][j] * (Sb[n] * O.Dkgrdpl * o.hg[j])
					o.P.dwgbdpg_n[i] += o.P.Mdl.Kgsat[i][j] * (Sb[n]*O.Dkgrdpg*o.hg[j] + kgr*dhgdpg_nj)
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
					o.P.Kll[m][n] -= coef * Gb[m][i] * o.P.dwlbdpl_n[i]
					o.P.Klg[m][n] -= coef * Gb[m][i] * o.P.dwlbdpg_n[i]
					o.P.Kgl[m][n] -= coef * Gb[m][i] * o.P.dwgbdpl_n[i]
					o.P.Kgg[m][n] -= coef * Gb[m][i] * o.P.dwgbdpg_n[i]
				}

				// extrapolation term
				if o.P.DoExtrap {
					o.P.dρldpl_ex[m][n] += o.P.Emat[m][idx] * O.Cpl * Sb[n]
					o.P.dρgdpg_ex[m][n] += o.P.Emat[m][idx] * O.Dpg * Sb[n]
				}
			}

			// extrapolation terms
			if o.P.DoExtrap {
				o.P.ρl_ex[n] += o.P.Emat[n][idx] * ρl
				o.P.ρg_ex[n] += o.P.Emat[n][idx] * ρg
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
			IpBmatrix(o.U.B, o.Ndim, u_nverts, G, radius, S, sol.Axisym)
			la.MatTrMulAdd3(o.U.K, coef, o.U.B, o.U.D, o.U.B) // K += coef * tr(B) * D * B
		} else {
			IpAddToKt(o.U.K, u_nverts, o.Ndim, coef, G, o.U.D)
		}
	}

	// contribution from natural boundary conditions
	if o.P.HasSeep {
		err = o.P.add_natbcs_to_jac(sol)
		if err != nil {
			return
		}
		err = o.add_natbcs_to_jac(sol)
		if err != nil {
			return
		}
	}

	// assemble K matrices into Kb
	o.P.assembleKs(Kb)
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
func (o *ElemUPP) Update(sol *Solution) (err error) {
	err = o.U.Update(sol)
	if err != nil {
		return
	}
	return o.P.Update(sol)
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// Ipoints returns the real coordinates of integration points [nip][ndim]
func (o *ElemUPP) Ipoints() (coords [][]float64) {
	coords = la.MatAlloc(len(o.U.IpsElem), o.Ndim)
	for idx, ip := range o.U.IpsElem {
		coords[idx] = o.U.Cell.Shp.IpRealCoords(o.U.X, ip)
	}
	return
}

// SetIniIvs sets initial ivs for given values in sol and ivs map
func (o *ElemUPP) SetIniIvs(sol *Solution, ivs map[string][]float64) (err error) {

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
func (o *ElemUPP) BackupIvs(aux bool) (err error) {
	err = o.U.BackupIvs(aux)
	if err != nil {
		return
	}
	return o.P.BackupIvs(aux)
}

// RestoreIvs restore internal variables from copies
func (o *ElemUPP) RestoreIvs(aux bool) (err error) {
	err = o.U.RestoreIvs(aux)
	if err != nil {
		return
	}
	return o.P.RestoreIvs(aux)
}

// Ureset fixes internal variables after u (displacements) have been zeroed
func (o *ElemUPP) Ureset(sol *Solution) (err error) {
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
func (o *ElemUPP) Encode(enc Encoder) (err error) {
	err = o.U.Encode(enc)
	if err != nil {
		return
	}
	return o.P.Encode(enc)
}

// Decode decodes internal variables
func (o *ElemUPP) Decode(dec Decoder) (err error) {
	err = o.U.Decode(dec)
	if err != nil {
		return
	}
	return o.P.Decode(dec)
}

// OutIpsData returns data from all integration points for output
func (o *ElemUPP) OutIpsData() (data []*OutIpData) {
	flowL := LiqFlowKeys(o.Ndim)
	flowG := GasFlowKeys(o.Ndim)
	sigs := StressKeys(o.Ndim)
	for idx, ip := range o.U.IpsElem {
		r := o.P.States[idx]
		s := o.U.States[idx]
		x := o.U.Cell.Shp.IpRealCoords(o.U.X, ip)
		calc := func(sol *Solution) (vals map[string]float64) {
			err := o.ipvars(idx, sol)
			if err != nil {
				return
			}
			ns := (1.0 - o.divus) * o.P.States[idx].A_ns0
			sl := r.A_sl
			sg := 1.0 - sl
			ρL := r.A_ρL
			ρG := r.A_ρG
			klr := o.P.Mdl.Cnd.Klr(sl)
			kgr := o.P.Mdl.Cnd.Kgr(sg)
			vals = map[string]float64{
				"sl": sl,
				"sg": sg,
				"pl": o.P.pl,
				"pg": o.P.pg,
				"pc": o.P.pg - o.P.pl,
				"nf": 1.0 - ns,
			}
			for i := 0; i < o.Ndim; i++ {
				for j := 0; j < o.Ndim; j++ {
					vals[flowL[i]] += klr * o.P.Mdl.Klsat[i][j] * o.hl[j] / ρL
					vals[flowG[i]] += kgr * o.P.Mdl.Kgsat[i][j] * o.hg[j] / ρG
				}
			}
			for i, _ := range sigs {
				vals[sigs[i]] = s.Sig[i]
			}
			return
		}
		data = append(data, &OutIpData{o.Id(), x, calc})
	}
	return
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// ipvars computes current values @ integration points. idx == index of integration point
func (o *ElemUPP) ipvars(idx int, sol *Solution) (err error) {

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
	o.P.compute_gvec(sol.T)

	// clear gpl and recover u-variables @ ip
	o.divus = 0
	for i := 0; i < o.Ndim; i++ {
		o.P.gpl[i], o.P.gpg[i], o.U.us[i] = 0, 0, 0
		for m := 0; m < o.U.Cell.Shp.Nverts; m++ {
			r := o.U.Umap[i+m*o.Ndim]
			o.U.us[i] += o.U.Cell.Shp.S[m] * sol.Y[r]
			o.divus += o.U.Cell.Shp.G[m][i] * sol.Y[r]
		}
	}

	// recover p-variables @ ip
	o.P.pl, o.P.pg = 0, 0
	for m := 0; m < o.P.Cell.Shp.Nverts; m++ {
		rl := o.P.Plmap[m]
		rg := o.P.Pgmap[m]
		o.P.pl += o.P.Cell.Shp.S[m] * sol.Y[rl]
		o.P.pg += o.P.Cell.Shp.S[m] * sol.Y[rg]
		for i := 0; i < o.Ndim; i++ {
			o.P.gpl[i] += o.P.Cell.Shp.G[m][i] * sol.Y[rl]
			o.P.gpg[i] += o.P.Cell.Shp.G[m][i] * sol.Y[rg]
		}
	}

	// compute b, hl and hg
	α1 := sol.DynCfs.α1
	for i := 0; i < o.Ndim; i++ {
		o.b[i] = α1*o.U.us[i] - o.U.ζs[idx][i] - o.P.g[i]
		o.hl[i] = -ρL*o.b[i] - o.P.gpl[i]
		o.hg[i] = -ρG*o.b[i] - o.P.gpg[i]
	}
	return
}

// add_natbcs_to_jac adds contribution from natural boundary conditions to Jacobian
func (o *ElemUPP) add_natbcs_to_jac(sol *Solution) (err error) {

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
				_, pl, fl = o.P.fipvars(iface, sol)
				plmax = o.P.Plmax[idx][jdx] - shift
				if plmax < 0 {
					plmax = 0
				}

				// compute derivatives
				g = pl - plmax // Eq. (24)
				rmp = o.P.ramp(fl + o.P.κ*g)
				for i, m := range o.P.Cell.Shp.FaceLocalVerts[iface] {
					for n := 0; n < u_nverts; n++ {
						for j := 0; j < o.Ndim; j++ {
							c := j + n*o.Ndim
							for l, r := range o.P.Cell.Shp.FaceLocalVerts[iface] {
								o.Klu[m][c] += coef * Sf[i] * Sf[l] * o.dρldus_ex[r][c] * rmp
							}
						}
					}
				}
			}
		}
	}
	return
}
