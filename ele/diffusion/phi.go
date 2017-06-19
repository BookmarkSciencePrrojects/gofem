// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package diffusion

import (
	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/shp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun/dbf"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// Phi implementes a general element to solve the following equation
//     dφ       ∂φ
//     -- + v . -- = s(x)
//     dt       ∂x
// Notes: v is a constant vector
type Phi struct {

	// basic data
	Cell *inp.Cell   // the cell structure
	X    [][]float64 // [ndim][nnode] matrix of nodal coordinates
	Nu   int         // total number of unknowns == number of vertices
	Ndim int         // space dimension

	// integration points
	IpsElem []shp.Ipoint // [nip] integration points of element

	// local starred variables
	Psi []float64 // [nip] ψ* = β1.φ + β2.dφdt

	// scratchpad. computed @ each ip
	K [][]float64 // [nu][nu] consistent tangent matrix

	// problem variables
	Umap []int // assembly map (location array/element equations)
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	ele.SetInfoFunc("phi", func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *ele.Info {

		// new info
		var info ele.Info

		nverts := cell.Shp.Nverts
		ykeys := []string{"h"}

		info.Dofs = make([][]string, nverts)
		for m := 0; m < nverts; m++ {
			info.Dofs[m] = ykeys
		}

		info.T1vars = ykeys

		// return information
		return &info
	})

	// element allocator
	ele.SetAllocator("phi", func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) ele.Element {

		// basic data
		var o Phi
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
		o.Psi = make([]float64, nip)

		// scratchpad. computed @ each ip
		o.K = la.MatAlloc(o.Nu, o.Nu)

		// return new element
		return &o
	})
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o *Phi) Id() int { return o.Cell.Id }

// SetEqs set equations
func (o *Phi) SetEqs(eqs [][]int, mixedform_eqs []int) (err error) {
	o.Umap = make([]int, o.Nu)
	for m := 0; m < o.Nu; m++ {
		o.Umap[m] = eqs[m][0]
	}
	return
}

// SetEleConds set element conditions
func (o *Phi) SetEleConds(key string, f dbf.T, extra string) (err error) {
	return
}

// InterpStarVars interpolate star variables to integration points
func (o *Phi) InterpStarVars(sol *ele.Solution) (err error) {

	// for each integration point
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		err = o.Cell.Shp.CalcAtIp(o.X, ip, false)
		if err != nil {
			return
		}

		//interpolate starred variables
		o.Psi[idx] = 0
		for m := 0; m < o.Nu; m++ {
			o.Psi[idx] += o.Cell.Shp.S[m] * sol.Psi[o.Umap[m]]
		}
	}
	return
}

// AddToRhs adds -R to global residual vector fb
func (o *Phi) AddToRhs(fb []float64, sol *ele.Solution) (err error) {

	// auxiliary
	β1 := sol.DynCfs.GetBet1()
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
func (o *Phi) AddToKb(Kb *la.Triplet, sol *ele.Solution, firstIt bool) (err error) {

	// auxiliary
	β1 := sol.DynCfs.GetBet1()
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
func (o *Phi) Encode(enc utl.Encoder) (err error) {
	return
}

// Decode decodes internal variables
func (o *Phi) Decode(dec utl.Decoder) (err error) {
	return
}

// OutIpCoords returns the coordinates of integration points
func (o *Phi) OutIpCoords() (C [][]float64) {
	C = make([][]float64, len(o.IpsElem))
	for idx, ip := range o.IpsElem {
		C[idx] = o.Cell.Shp.IpRealCoords(o.X, ip)
	}
	return
}

// OutIpKeys returns the integration points' keys
func (o *Phi) OutIpKeys() []string {
	if o.Ndim == 3 {
		return []string{"Gphix", "Gphiy", "Gphiz"}
	}
	return []string{"Gphix", "Gphiy"}
}

// OutIpVals returns the integration points' values corresponding to keys
func (o *Phi) OutIpVals(M *ele.IpsMap, sol *ele.Solution) {
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
