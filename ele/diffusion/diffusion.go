// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package diffusion implements elements for diffusion problems
package diffusion

import (
	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mdl/diffusion"
	"github.com/cpmech/gofem/shp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// Diffusion implements an element for solving the diffusion equation expressed as
//
//     du                                      du
//   ρ ── + div w = s      with      w = -k(u) ──
//     dt                                      dx
//
type Diffusion struct {

	// basic data
	Cell *inp.Cell     // the cell structure
	X    [][]float64   // matrix of nodal coordinates [ndim][nnode]
	Ndim int           // space dimension
	Umap []int         // assembly map (location array/element equations)
	Mdl  *diffusion.M1 // model
	Sfun fun.Func      // s(x) function

	// integration points
	IpsElem []shp.Ipoint // integration points of element
	IpsFace []shp.Ipoint // integration points corresponding to faces

	// natural boundary conditions
	NatBcs []*ele.NaturalBc // natural boundary conditions

	// scratchpad
	Ustar []float64   // [nip] ustar* = β1 u + β2 dudt
	Xip   []float64   // real coordinates of ip
	Uval  float64     // u(t,x) scalar field @ ip
	Gradu []float64   // [ndim] ∇u(t,x): gradient of u @ ip
	Wvec  []float64   // [ndim] w(t,x) vector @ ip
	Tmp   []float64   // auxiliary vector
	K     [][]float64 // Jacobian matrix
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	ele.SetInfoFunc("diffusion", func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *ele.Info {

		// new info
		var info ele.Info
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
	})

	// element allocator
	ele.SetAllocator("diffusion", func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) ele.Element {

		// basic data
		var o Diffusion
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
		o.Mdl = mat.Dif.(*diffusion.M1)
		if sim.Data.Steady {
			o.Mdl.Rho = 0
		}

		// set natural boundary conditions
		for _, fc := range cell.FaceBcs {
			o.NatBcs = append(o.NatBcs, &ele.NaturalBc{fc.Cond, fc.FaceId, fc.Func, fc.Extra})
		}

		// scratchpad
		nverts := cell.Shp.Nverts
		o.Ustar = make([]float64, nip)
		o.Xip = make([]float64, o.Ndim)
		o.Gradu = make([]float64, o.Ndim)
		o.Wvec = make([]float64, o.Ndim)
		o.Tmp = make([]float64, o.Ndim)
		o.K = la.MatAlloc(nverts, nverts)
		return &o
	})
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o *Diffusion) Id() int { return o.Cell.Id }

// SetEqs sets equations
func (o *Diffusion) SetEqs(eqs [][]int, mixedform_eqs []int) (err error) {
	nverts := o.Cell.Shp.Nverts
	o.Umap = make([]int, nverts)
	for m := 0; m < nverts; m++ {
		o.Umap[m] = eqs[m][0]
	}
	return
}

// SetEleConds sets element conditions
func (o *Diffusion) SetEleConds(key string, f fun.Func, extra string) (err error) {
	o.Sfun = nil
	if key == "s" {
		o.Sfun = f
	}
	return
}

// InterpStarVars interpolates star variables to integration points
func (o *Diffusion) InterpStarVars(sol *ele.Solution) (err error) {

	// for each integration point
	nverts := o.Cell.Shp.Nverts
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}

		// interpolate starred variables
		o.Ustar[idx] = 0
		for m := 0; m < nverts; m++ {
			o.Ustar[idx] += o.Cell.Shp.S[m] * sol.Psi[o.Umap[m]]
		}
	}
	return
}

// AddToRhs adds -R to global residual vector fb
func (o *Diffusion) AddToRhs(fb []float64, sol *ele.Solution) (err error) {

	// for each integration point
	ρ := o.Mdl.Rho
	β1 := sol.DynCfs.GetBet1()
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
		dudt = β1*o.Uval - o.Ustar[idx]
		kval = o.Mdl.Kval(o.Uval)
		if o.Sfun != nil {
			sval = o.Sfun.F(sol.T, o.Xip)
		}

		// compute Wvec
		for i := 0; i < o.Ndim; i++ {
			o.Wvec[i] = 0
			for j := 0; j < o.Ndim; j++ {
				o.Wvec[i] -= kval * o.Mdl.Kcte[i][j] * o.Gradu[j]
			}
		}

		// add negative of residual term to fb
		for m := 0; m < nverts; m++ {
			r := o.Umap[m]
			fb[r] -= coef * S[m] * (ρ*dudt - sval) // - ftrs + fext/sval
			for i := 0; i < o.Ndim; i++ {
				fb[r] += coef * G[m][i] * o.Wvec[i] // + fint
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
func (o *Diffusion) AddToKb(Kb *la.Triplet, sol *ele.Solution, firstIt bool) (err error) {

	// clear matrices
	la.MatFill(o.K, 0)

	// for each integration point
	ρ := o.Mdl.Rho
	β1 := sol.DynCfs.GetBet1()
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
		kval = o.Mdl.Kval(o.Uval)
		dkdu = o.Mdl.DkDu(o.Uval)

		// K := dR/du
		for n := 0; n < nverts; n++ {
			for j := 0; j < o.Ndim; j++ {
				o.Tmp[j] = S[n]*dkdu*o.Gradu[j] + kval*G[n][j]
			}
			for m := 0; m < nverts; m++ {
				o.K[m][n] += coef * S[m] * S[n] * β1 * ρ
				for i := 0; i < o.Ndim; i++ {
					for j := 0; j < o.Ndim; j++ {
						o.K[m][n] += coef * G[m][i] * o.Mdl.Kcte[i][j] * o.Tmp[j]
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
func (o *Diffusion) Update(sol *ele.Solution) (err error) {
	return
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// SetIniIvs sets initial ivs for given values in sol and ivs map
func (o *Diffusion) SetIniIvs(sol *ele.Solution, ignored map[string][]float64) (err error) {
	return
}

// BackupIvs creates copy of internal variables
func (o *Diffusion) BackupIvs(aux bool) (err error) {
	return
}

// RestoreIvs restores internal variables from copies
func (o *Diffusion) RestoreIvs(aux bool) (err error) {
	return
}

// Ureset fixes internal variables after u (displacements) have been zeroed
func (o *Diffusion) Ureset(sol *ele.Solution) (err error) {
	return
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o *Diffusion) Encode(enc utl.Encoder) (err error) {
	return
}

// Decode decodes internal variables
func (o *Diffusion) Decode(dec utl.Decoder) (err error) {
	return
}

// OutIpCoords returns the coordinates of integration points
func (o *Diffusion) OutIpCoords() (C [][]float64) {
	C = make([][]float64, len(o.IpsElem))
	for idx, ip := range o.IpsElem {
		C[idx] = o.Cell.Shp.IpRealCoords(o.X, ip)
	}
	return
}

// OutIpKeys returns the integration points' keys
func (o *Diffusion) OutIpKeys() []string {
	return nil
}

// OutIpVals returns the integration points' values corresponding to keys
func (o *Diffusion) OutIpVals(M *ele.IpsMap, sol *ele.Solution) {
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// ipvars computes current values @ integration points. idx == index of integration point
func (o *Diffusion) ipvars(idx int, sol *ele.Solution) (err error) {

	// interpolation functions and gradients
	err = o.Cell.Shp.CalcAtIp(o.X, o.IpsElem[idx], true)
	if err != nil {
		return
	}

	// clear Uval and its gradient @ ip
	o.Uval = 0
	for i := 0; i < o.Ndim; i++ {
		o.Gradu[i] = 0
		o.Xip[i] = 0
	}

	// compute u and its gradient @ ip by means of interpolating from nodes
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		r := o.Umap[m]
		o.Uval += o.Cell.Shp.S[m] * sol.Y[r]
		for i := 0; i < o.Ndim; i++ {
			o.Gradu[i] += o.Cell.Shp.G[m][i] * sol.Y[r]
			o.Xip[i] += o.Cell.Shp.S[m] * o.X[i][m]
		}
	}
	return
}

// add_natbcs_to_rhs adds natural boundary conditions to rhs
func (o *Diffusion) add_natbcs_to_rhs(fb []float64, sol *ele.Solution) (err error) {

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
