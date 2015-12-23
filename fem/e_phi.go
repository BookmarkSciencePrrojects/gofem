// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/shp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
)

// ElemPhi implementes a general element to solve the following equation
//     dφ       ∂φ
//     -- + v . -- = s(x)
//     dt       ∂x
// Notes: v is a constant vector
type ElemPhi struct {

	// basic data
	Cell *inp.Cell   // the cell structure
	X    [][]float64 // [ndim][nnode] matrix of nodal coordinates
	Nu   int         // total number of unknowns == number of vertices
	Ndim int         // space dimension

	// integration points
	IpsElem []shp.Ipoint // [nip] integration points of element

	// local starred variables
	ψs []float64 // [nip] ψ* = β1.φ + β2.dφdt

	// scratchpad. computed @ each ip
	K [][]float64 // [nu][nu] consistent tangent matrix

	// problem variables
	Umap []int // assembly map (location array/element equations)
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	infogetters["phi"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *Info {

		// new info
		var info Info

		nverts := cell.Shp.Nverts
		ykeys := []string{"h"}

		info.Dofs = make([][]string, nverts)
		for m := 0; m < nverts; m++ {
			info.Dofs[m] = ykeys
		}

		info.T1vars = ykeys

		// return information
		return &info
	}

	// element allocator
	eallocators["phi"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) Elem {

		// basic data
		var o ElemPhi
		o.Cell = cell
		o.X = x
		o.Nu = o.Cell.Shp.Nverts
		o.Ndim = sim.Ndim

		// integration points
		var err error
		o.IpsElem, _, err = o.Cell.Shp.GetIps(edat.Nip, edat.Nipf)
		if err != nil {
			chk.Panic("cannot allocate integration points of phi-element with nip=%d and nipf=%d:\n%v", edat.Nip, edat.Nipf, err)
		}

		// local starred variables
		nip := len(o.IpsElem)
		o.ψs = make([]float64, nip)

		// scratchpad. computed @ each ip
		o.K = la.MatAlloc(o.Nu, o.Nu)

		// return new element
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o *ElemPhi) Id() int { return o.Cell.Id }

// SetEqs set equations
func (o *ElemPhi) SetEqs(eqs [][]int, mixedform_eqs []int) (err error) {
	o.Umap = make([]int, o.Nu)
	for m := 0; m < o.Nu; m++ {
		o.Umap[m] = eqs[m][0]
	}
	return
}

// SetEleConds set element conditions
func (o *ElemPhi) SetEleConds(key string, f fun.Func, extra string) (err error) {
	return
}

// InterpStarVars interpolate star variables to integration points
func (o *ElemPhi) InterpStarVars(sol *Solution) (err error) {

	// for each integration point
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, false)
		if err != nil {
			return
		}

		//interpolate starred variables
		o.ψs[idx] = 0
		for m := 0; m < o.Nu; m++ {
			o.ψs[idx] += o.Cell.Shp.S[m] * sol.Psi[o.Umap[m]]
		}
	}
	return
}

// AddToRhs adds -R to global residual vector fb
func (o *ElemPhi) AddToRhs(fb []float64, sol *Solution) (err error) {

	// auxiliary
	β1 := sol.DynCfs.β1
	nverts := o.Cell.Shp.Nverts

	v := []float64{1, 0, 0} // TODO: find a way to input the velocity

	// for each integration point
	for _, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}

		// auxiliary variables
		coef := o.Cell.Shp.J * ip[3]
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G

		// add to right hand side vector
		for m := 0; m < nverts; m++ {
			r := o.Umap[m] // row in the global vector
			for n := 0; n < nverts; n++ {
				fb[r] -= coef * S[m] * S[n] * β1 * sol.Y[o.Umap[n]]
				for j := 0; j < o.Ndim; j++ {
					fb[r] += coef * v[j] * G[m][j] * S[n]
				}
			}
		}
	}
	return
}

// AddToKb adds element K to global Jacobian matrix Kb
func (o *ElemPhi) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (err error) {

	// auxiliary
	β1 := sol.DynCfs.β1
	nverts := o.Cell.Shp.Nverts

	// zero K matrix
	la.MatFill(o.K, 0)

	// for each integration point
	for _, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}

		// auxiliary variables
		coef := o.Cell.Shp.J * ip[3]
		S := o.Cell.Shp.S

		// add to right hand side vector
		for m := 0; m < nverts; m++ {
			for n := 0; n < nverts; n++ {
				o.K[m][n] += coef * S[m] * S[n] * β1
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

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o *ElemPhi) Encode(enc Encoder) (err error) {
	return
}

// Decode decodes internal variables
func (o *ElemPhi) Decode(dec Decoder) (err error) {
	return
}

// OutIpCoords returns the coordinates of integration points
func (o *ElemPhi) OutIpCoords() (C [][]float64) {
	C = make([][]float64, len(o.IpsElem))
	for idx, ip := range o.IpsElem {
		C[idx] = o.Cell.Shp.IpRealCoords(o.X, ip)
	}
	return
}

// OutIpKeys returns the integration points' keys
func (o *ElemPhi) OutIpKeys() []string {
	if o.Ndim == 3 {
		return []string{"Gphix", "Gphiy", "Gphiz"}
	}
	return []string{"Gphix", "Gphiy"}
}

// OutIpVals returns the integration points' values corresponding to keys
func (o *ElemPhi) OutIpVals(M *IpsMap, sol *Solution) {
	keys := o.OutIpKeys()
	nip := len(o.IpsElem)
	for idx, ip := range o.IpsElem {
		err := o.Cell.Shp.CalcAtIp(o.X, ip, true)
		if err != nil {
			return
		}
		G := o.Cell.Shp.G
		for i := 0; i < o.Ndim; i++ {
			var Gphi_i float64
			for m := 0; m < o.Nu; m++ {
				Gphi_i += G[m][i] * sol.Y[o.Umap[m]]
			}
			M.Set(keys[i], idx, nip, Gphi_i)
		}
	}
}
