// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mdl/sld"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
)

// Rod represents a structural rod element (for only axial loads)
type Rod struct {

	// basic data
	Cell *inp.Cell   // the cell structure
	X    [][]float64 // matrix of nodal coordinates [ndim][nnode]
	Nu   int         // total number of unknowns == 2 * nsn
	Ndim int         // space dimension

	// variables for dynamics
	Gfcn fun.Func // gravity function

	// integration points
	IpsElem []shp.Ipoint // integration points of element

	// vectors and matrices
	K [][]float64 // element K matrix
	M [][]float64 // element M matrix

	// problem variables
	Umap []int // assembly map (location array/element equations)

	// material model and internal variables
	Mdl       sld.OneD
	States    []*sld.OnedState
	StatesBkp []*sld.OnedState
	StatesAux []*sld.OnedState

	// scratchpad. computed @ each ip
	grav []float64 // [ndim] gravity vector
	us   []float64 // [ndim] displacements @ ip
	fi   []float64 // [nu] internal forces
	ue   []float64 // local u vector
}

// register element
func init() {

	// information allocator
	infogetters["rod"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *Info {

		// new info
		var info Info

		// number of nodes in element
		nverts := cell.Shp.Nverts

		// solution variables
		ykeys := []string{"ux", "uy"}
		if sim.Ndim == 3 {
			ykeys = []string{"ux", "uy", "uz"}
		}
		info.Dofs = make([][]string, nverts)
		for m := 0; m < nverts; m++ {
			info.Dofs[m] = ykeys
		}

		// maps
		info.Y2F = map[string]string{"ux": "fx", "uy": "fy", "uz": "fz"}

		// t1 and t2 variables
		info.T2vars = ykeys
		return &info
	}

	// element allocator
	eallocators["rod"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) Elem {

		// basic data
		var o Rod
		o.Cell = cell
		o.X = x
		o.Ndim = sim.Ndim
		o.Nu = o.Ndim * o.Cell.Shp.Nverts

		// model
		mat := sim.MatModels.Get(edat.Mat)
		if mat == nil {
			chk.Panic("cannot find material %q for Rod {tag=%d, id=%d}\n", edat.Mat, cell.Tag, cell.Id)
		}
		o.Mdl = mat.Sld.(sld.OneD)

		// integration points
		var err error
		o.IpsElem, _, err = o.Cell.Shp.GetIps(edat.Nip, 0)
		if err != nil {
			chk.Panic("cannot get integration points for rod element {tag=%d id=%d material=%q} with nip=%d", cell.Tag, cell.Id, edat.Mat, edat.Nip)
		}

		// scratchpad. computed @ each ip
		o.K = la.MatAlloc(o.Nu, o.Nu)
		o.M = la.MatAlloc(o.Nu, o.Nu)
		o.ue = make([]float64, o.Nu)

		// scratchpad. computed @ each ip
		o.grav = make([]float64, o.Ndim)
		o.us = make([]float64, o.Ndim)
		o.fi = make([]float64, o.Nu)

		// return new element
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o *Rod) Id() int { return o.Cell.Id }

// SetEqs set equations
func (o *Rod) SetEqs(eqs [][]int, mixedform_eqs []int) (err error) {
	o.Umap = make([]int, o.Nu)
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		for i := 0; i < o.Ndim; i++ {
			r := i + m*o.Ndim
			o.Umap[r] = eqs[m][i]
		}
	}
	return
}

// InterpStarVars interpolates star variables to integration points
func (o *Rod) InterpStarVars(sol *Solution) (err error) {
	return
}

// SetEleConds set element conditions
func (o *Rod) SetEleConds(key string, f fun.Func, extra string) (err error) {
	if key == "g" {
		o.Gfcn = f
	}
	return
}

// AddToRhs adds -R to global residual vector fb
func (o *Rod) AddToRhs(fb []float64, sol *Solution) (err error) {

	// for each integration point
	A := o.Mdl.GetA()
	nverts := o.Cell.Shp.Nverts
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		err = o.ipvars(idx, sol)
		if err != nil {
			return
		}

		// auxiliary
		coef := ip[3]
		Jvec := o.Cell.Shp.Jvec3d
		G := o.Cell.Shp.Gvec
		σ := o.States[idx].Sig

		// update fb with internal forces
		for m := 0; m < nverts; m++ {
			for i := 0; i < o.Ndim; i++ {
				r := o.Umap[i+m*o.Ndim]
				fb[r] -= coef * A * σ * G[m] * Jvec[i] // -fi
			}
		}
	}
	return
}

// AddToKb adds element K to global Jacobian matrix Kb
func (o *Rod) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (err error) {

	// zero K matrix
	la.MatFill(o.K, 0)
	la.MatFill(o.M, 0) // TODO: implement mass matrix

	// for each integration point
	A := o.Mdl.GetA()
	var E float64
	nverts := o.Cell.Shp.Nverts
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		err = o.ipvars(idx, sol)
		if err != nil {
			return
		}

		// auxiliary
		coef := ip[3]
		Jvec := o.Cell.Shp.Jvec3d
		G := o.Cell.Shp.Gvec
		J := o.Cell.Shp.J

		// add contribution to consistent tangent matrix
		for m := 0; m < nverts; m++ {
			for n := 0; n < nverts; n++ {
				for i := 0; i < o.Ndim; i++ {
					for j := 0; j < o.Ndim; j++ {
						r := i + m*o.Ndim
						c := j + n*o.Ndim
						E, _, err = o.Mdl.CalcD(o.States[idx], firstIt)
						if err != nil {
							return
						}
						o.K[r][c] += coef * A * E * G[m] * G[n] * Jvec[i] * Jvec[j] / J
					}
				}
			}
		}
	}

	// add K to sparse matrix Kb
	for i, I := range o.Umap {
		for j, J := range o.Umap {
			Kb.Put(I, J, o.K[i][j])
		}
	}
	return
}

// Update perform (tangent) update
func (o *Rod) Update(sol *Solution) (err error) {

	// for each integration point
	nverts := o.Cell.Shp.Nverts
	for idx, _ := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		err = o.ipvars(idx, sol)
		if err != nil {
			return
		}

		// auxiliary
		Jvec := o.Cell.Shp.Jvec3d
		G := o.Cell.Shp.Gvec
		J := o.Cell.Shp.J

		// compute strains
		Δε := 0.0
		for m := 0; m < nverts; m++ {
			for i := 0; i < o.Ndim; i++ {
				r := o.Umap[i+m*o.Ndim]
				Δε += G[m] * Jvec[i] * sol.ΔY[r] / J
			}
		}

		// call model update => update stresses
		err = o.Mdl.Update(o.States[idx], 0.0, Δε, 0)
		if err != nil {
			return
		}
	}
	return
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// Ipoints returns the real coordinates of integration points [nip][ndim]
func (o *Rod) Ipoints() (coords [][]float64) {
	coords = la.MatAlloc(len(o.IpsElem), o.Ndim)
	for idx, ip := range o.IpsElem {
		coords[idx] = o.Cell.Shp.IpRealCoords(o.X, ip)
	}
	return
}

// SetIniIvs sets initial ivs for given values in sol and ivs map
func (o *Rod) SetIniIvs(sol *Solution, ivs map[string][]float64) (err error) {

	// allocate slices of states
	nip := len(o.IpsElem)
	o.States = make([]*sld.OnedState, nip)
	o.StatesBkp = make([]*sld.OnedState, nip)
	o.StatesAux = make([]*sld.OnedState, nip)

	// for each integration point
	for i := 0; i < nip; i++ {
		o.States[i], _ = o.Mdl.InitIntVars1D()
		o.StatesBkp[i] = o.States[i].GetCopy()
		o.StatesAux[i] = o.States[i].GetCopy()
	}

	// initial stresses
	if _, ok := ivs["sig"]; ok {
		for i := 0; i < nip; i++ {
			o.States[i].Sig = ivs["sig"][i]
			o.StatesBkp[i].Sig = o.States[i].Sig
		}
	}
	return
}

// SetIvs set secondary variables; e.g. during initialisation via files
func (o *Rod) SetIvs(zvars map[string][]float64) (err error) {
	return
}

// BackupIvs create copy of internal variables
func (o *Rod) BackupIvs(aux bool) (err error) {
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
func (o *Rod) RestoreIvs(aux bool) (err error) {
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
func (o *Rod) Ureset(sol *Solution) (err error) {
	return
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o *Rod) Encode(enc Encoder) (err error) {
	return enc.Encode(o.States)
}

// Decode decodes internal variables
func (o *Rod) Decode(dec Decoder) (err error) {
	err = dec.Decode(&o.States)
	if err != nil {
		return
	}
	return o.BackupIvs(false)
}

// OutIpCoords returns the coordinates of integration points
func (o *Rod) OutIpCoords() (C [][]float64) {
	C = make([][]float64, len(o.IpsElem))
	for idx, ip := range o.IpsElem {
		C[idx] = o.Cell.Shp.IpRealCoords(o.X, ip)
	}
	return
}

// OutIpKeys returns the integration points' keys
func (o *Rod) OutIpKeys() []string {
	return []string{"sig"}
}

// OutIpVals returns the integration points' values corresponding to keys
func (o *Rod) OutIpVals(M *IpsMap, sol *Solution) {
	nip := len(o.IpsElem)
	for idx, _ := range o.IpsElem {
		M.Set("sig", idx, nip, o.States[idx].Sig)
	}
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// ipvars computes current values @ integration points. idx == index of integration point
func (o *Rod) ipvars(idx int, sol *Solution) (err error) {

	// interpolation functions and gradients
	err = o.Cell.Shp.CalcAtIp(o.X, o.IpsElem[idx], true)
	if err != nil {
		return
	}

	// skip if steady (this must be after CalcAtIp, because callers will need S and G)
	if sol.Steady {
		return
	}

	// clear variables
	for i := 0; i < o.Ndim; i++ {
		o.us[i] = 0
	}

	// recover u-variables @ ip
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		for i := 0; i < o.Ndim; i++ {
			r := o.Umap[i+m*o.Ndim]
			o.us[i] += o.Cell.Shp.S[m] * sol.Y[r]
		}
	}
	return
}
