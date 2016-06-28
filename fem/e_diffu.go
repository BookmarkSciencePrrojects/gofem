// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mdl/gen"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
)

// ElemDiffu implements an element for solving the diffusion equation expressed as
//
//     du                                      du
//   ρ ── + div w = s      with      w = -k(u) ──
//     dt                                      dx
//
type ElemDiffu struct {

	// basic data
	Cell *inp.Cell    // the cell structure
	X    [][]float64  // matrix of nodal coordinates [ndim][nnode]
	Ndim int          // space dimension
	Umap []int        // assembly map (location array/element equations)
	Mdl  *gen.Diffu01 // model
	Sfun fun.Func     // s(x) function

	// integration points
	IpsElem []shp.Ipoint // integration points of element
	IpsFace []shp.Ipoint // integration points corresponding to faces

	// natural boundary conditions
	NatBcs []*NaturalBc // natural boundary conditions

	// scratchpad
	ustar []float64   // [nip] ustar* = β1 u + β2 dudt
	xip   []float64   // real coordinates of ip
	uval  float64     // u(t,x) scalar field @ ip
	gradu []float64   // [ndim] ∇u(t,x): gradient of u @ ip
	wvec  []float64   // [ndim] w(t,x) vector @ ip
	tmp   []float64   // auxiliary vector
	K     [][]float64 // Jacobian matrix
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	infogetters["diffu"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *Info {

		// new info
		var info Info
		nverts := cell.Shp.Nverts

		// solution variables
		ykeys := []string{"u"}
		info.Dofs = make([][]string, nverts)
		for m := 0; m < nverts; m++ {
			info.Dofs[m] = ykeys
		}

		// Y2F map and t1 and t2 variables
		info.Y2F = map[string]string{"u": "q"}
		info.T1vars = ykeys
		return &info
	}

	// element allocator
	eallocators["diffu"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) Elem {

		// basic data
		var o ElemDiffu
		o.Cell = cell
		o.X = x
		o.Ndim = sim.Ndim

		// integration points
		var err error
		o.IpsElem, o.IpsFace, err = o.Cell.Shp.GetIps(edat.Nip, edat.Nipf)
		if err != nil {
			chk.Panic("cannot allocate integration points of diffusion element with nip=%d and nipf=%d:\n%v", edat.Nip, edat.Nipf, err)
		}
		nip := len(o.IpsElem)

		// model
		mat := sim.MatModels.Get(edat.Mat)
		if mat == nil {
			chk.Panic("cannot get model for diffusion element {tag=%d id=%d material=%q}:\n%v", cell.Tag, cell.Id, edat.Mat, err)
		}
		o.Mdl = mat.Gen.(*gen.Diffu01)
		if sim.Data.Steady {
			o.Mdl.Rho = 0
		}

		// set natural boundary conditions
		for _, fc := range cell.FaceBcs {
			o.NatBcs = append(o.NatBcs, &NaturalBc{fc.Cond, fc.FaceId, fc.Func, fc.Extra})
		}

		// scratchpad
		nverts := cell.Shp.Nverts
		o.ustar = make([]float64, nip)
		o.xip = make([]float64, o.Ndim)
		o.gradu = make([]float64, o.Ndim)
		o.wvec = make([]float64, o.Ndim)
		o.tmp = make([]float64, o.Ndim)
		o.K = la.MatAlloc(nverts, nverts)
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o *ElemDiffu) Id() int { return o.Cell.Id }

// SetEqs sets equations
func (o *ElemDiffu) SetEqs(eqs [][]int, mixedform_eqs []int) (err error) {
	nverts := o.Cell.Shp.Nverts
	o.Umap = make([]int, nverts)
	for m := 0; m < nverts; m++ {
		o.Umap[m] = eqs[m][0]
	}
	return
}

// SetEleConds sets element conditions
func (o *ElemDiffu) SetEleConds(key string, f fun.Func, extra string) (err error) {
	o.Sfun = nil
	if key == "s" {
		o.Sfun = f
	}
	return
}

// InterpStarVars interpolates star variables to integration points
func (o *ElemDiffu) InterpStarVars(sol *Solution) (err error) {

	// for each integration point
	nverts := o.Cell.Shp.Nverts
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}

		// interpolate starred variables
		o.ustar[idx] = 0
		for m := 0; m < nverts; m++ {
			o.ustar[idx] += o.Cell.Shp.S[m] * sol.Psi[o.Umap[m]]
		}
	}
	return
}

// AddToRhs adds -R to global residual vector fb
func (o *ElemDiffu) AddToRhs(fb []float64, sol *Solution) (err error) {

	// for each integration point
	ρ := o.Mdl.Rho
	β1 := sol.DynCfs.β1
	nverts := o.Cell.Shp.Nverts
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
		dudt = β1*o.uval - o.ustar[idx]
		kval = o.Mdl.Kval(o.uval)
		if o.Sfun != nil {
			sval = o.Sfun.F(sol.T, o.xip)
		}

		// compute wvec
		for i := 0; i < o.Ndim; i++ {
			o.wvec[i] = 0
			for j := 0; j < o.Ndim; j++ {
				o.wvec[i] -= kval * o.Mdl.Kcte[i][j] * o.gradu[j]
			}
		}

		// add negative of residual term to fb
		for m := 0; m < nverts; m++ {
			r := o.Umap[m]
			fb[r] -= coef * S[m] * (ρ*dudt - sval) // - ftrs + fext/sval
			for i := 0; i < o.Ndim; i++ {
				fb[r] += coef * G[m][i] * o.wvec[i] // + fint
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
func (o *ElemDiffu) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (err error) {

	// clear matrices
	la.MatFill(o.K, 0)

	// for each integration point
	ρ := o.Mdl.Rho
	β1 := sol.DynCfs.β1
	nverts := o.Cell.Shp.Nverts
	var coef, kval, dkdu float64
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		err = o.ipvars(idx, sol)
		if err != nil {
			return
		}
		coef = o.Cell.Shp.J * ip[3]
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G
		kval = o.Mdl.Kval(o.uval)
		dkdu = o.Mdl.DkDu(o.uval)

		// K := dR/du
		for n := 0; n < nverts; n++ {
			for j := 0; j < o.Ndim; j++ {
				o.tmp[j] = S[n]*dkdu*o.gradu[j] + kval*G[n][j]
			}
			for m := 0; m < nverts; m++ {
				o.K[m][n] += coef * S[m] * S[n] * β1 * ρ
				for i := 0; i < o.Ndim; i++ {
					for j := 0; j < o.Ndim; j++ {
						o.K[m][n] += coef * G[m][i] * o.Mdl.Kcte[i][j] * o.tmp[j]
					}
				}
			}
		}
	}

	// add to sparse matrix Kb
	for i, I := range o.Umap {
		for j, J := range o.Umap {
			Kb.Put(I, J, o.K[i][j])
		}
	}
	return
}

// Update performs (tangent) update
func (o *ElemDiffu) Update(sol *Solution) (err error) {
	return
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// SetIniIvs sets initial ivs for given values in sol and ivs map
func (o *ElemDiffu) SetIniIvs(sol *Solution, ignored map[string][]float64) (err error) {
	return
}

// BackupIvs creates copy of internal variables
func (o *ElemDiffu) BackupIvs(aux bool) (err error) {
	return
}

// RestoreIvs restores internal variables from copies
func (o *ElemDiffu) RestoreIvs(aux bool) (err error) {
	return
}

// Ureset fixes internal variables after u (displacements) have been zeroed
func (o *ElemDiffu) Ureset(sol *Solution) (err error) {
	return
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o *ElemDiffu) Encode(enc Encoder) (err error) {
	return
}

// Decode decodes internal variables
func (o *ElemDiffu) Decode(dec Decoder) (err error) {
	return
}

// OutIpCoords returns the coordinates of integration points
func (o *ElemDiffu) OutIpCoords() (C [][]float64) {
	C = make([][]float64, len(o.IpsElem))
	for idx, ip := range o.IpsElem {
		C[idx] = o.Cell.Shp.IpRealCoords(o.X, ip)
	}
	return
}

// OutIpKeys returns the integration points' keys
func (o *ElemDiffu) OutIpKeys() []string {
	return nil
}

// OutIpVals returns the integration points' values corresponding to keys
func (o *ElemDiffu) OutIpVals(M *IpsMap, sol *Solution) {
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// ipvars computes current values @ integration points. idx == index of integration point
func (o *ElemDiffu) ipvars(idx int, sol *Solution) (err error) {

	// interpolation functions and gradients
	err = o.Cell.Shp.CalcAtIp(o.X, o.IpsElem[idx], true)
	if err != nil {
		return
	}

	// clear uval and its gradient @ ip
	o.uval = 0
	for i := 0; i < o.Ndim; i++ {
		o.gradu[i] = 0
		o.xip[i] = 0
	}

	// compute u and its gradient @ ip by means of interpolating from nodes
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		r := o.Umap[m]
		o.uval += o.Cell.Shp.S[m] * sol.Y[r]
		for i := 0; i < o.Ndim; i++ {
			o.gradu[i] += o.Cell.Shp.G[m][i] * sol.Y[r]
			o.xip[i] += o.Cell.Shp.S[m] * o.X[i][m]
		}
	}
	return
}

// add_natbcs_to_rhs adds natural boundary conditions to rhs
func (o *ElemDiffu) add_natbcs_to_rhs(fb []float64, sol *Solution) (err error) {

	// compute surface integral
	var qb float64
	for _, nbc := range o.NatBcs {

		// specified flux
		qb = nbc.Fcn.F(sol.T, nil)

		// loop over ips of face
		for _, ipf := range o.IpsFace {

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
			case "qb":
				for i, m := range o.Cell.Shp.FaceLocalVerts[iface] {
					fb[o.Umap[m]] -= coef * qb * Sf[i]
				}
			}
		}
	}
	return
}
