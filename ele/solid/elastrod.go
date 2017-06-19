// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package solid

import (
	"math"

	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mdl/solid"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun/dbf"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// ElastRod represents a structural rod element (for axial loads only) with 2 nodes only and
// simply implemented with constant stiffness matrix; i.e. no numerical integration is needed
type ElastRod struct {

	// basic data
	Cell *inp.Cell   // the cell structure
	X    [][]float64 // matrix of nodal coordinates [ndim][nnode]
	Nu   int         // total number of unknowns == 2 * nsn
	Ndim int         // space dimension

	// parameters and properties
	Mdl *solid.OnedLinElast // material model with: E, G, A, I22, I11, Jtt and Rho
	L   float64             // length of rod

	// variables for dynamics
	Gfcn dbf.T // gravity function

	// vectors and matrices
	T [][]float64 // [ndim][nu] transformation matrix: system aligned to rod => element system
	K [][]float64 // [nu][nu] element K matrix
	M [][]float64 // [nu][nu] element M matrix

	// problem variables
	Umap []int // assembly map (location array/element equations)

	// scratchpad. computed @ each ip
	ua []float64 // [2] local axial displacements
}

// register element
func init() {

	// information allocator
	ele.SetInfoFunc("elastrod", func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *ele.Info {

		// new info
		var info ele.Info

		// solution variables
		ykeys := []string{"ux", "uy"}
		if sim.Ndim == 3 {
			ykeys = []string{"ux", "uy", "uz"}
		}
		info.Dofs = make([][]string, 2)
		for m := 0; m < 2; m++ {
			info.Dofs[m] = ykeys
		}

		// maps
		info.Y2F = map[string]string{"ux": "fx", "uy": "fy", "uz": "fz"}

		// t1 and t2 variables
		info.T2vars = ykeys
		return &info
	})

	// element allocator
	ele.SetAllocator("elastrod", func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) ele.Element {

		// check
		ndim := len(x)
		if ndim == 3 {
			chk.Panic("elastrod is not implemented for 3D yet")
		}

		// basic data
		var o ElastRod
		o.Cell = cell
		o.X = x
		o.Ndim = sim.Ndim
		o.Nu = o.Ndim * 2

		// parameters
		mat := sim.MatModels.Get(edat.Mat)
		if mat == nil {
			chk.Panic("cannot get materials data for elastic rod element {tag=%d id=%d material=%q}", cell.Tag, cell.Id, edat.Mat)
		}
		o.Mdl = mat.Sld.(*solid.OnedLinElast)

		// vectors and matrices
		o.K = la.MatAlloc(o.Nu, o.Nu)
		o.M = la.MatAlloc(o.Nu, o.Nu)
		o.ua = make([]float64, 2)

		// geometry
		x0 := o.X[0][0]
		y0 := o.X[1][0]
		x1 := o.X[0][1]
		y1 := o.X[1][1]
		dx := x1 - x0
		dy := y1 - y0
		o.L = math.Sqrt(dx*dx + dy*dy)

		// global-to-local transformation matrix
		c := dx / o.L
		s := dy / o.L
		o.T = [][]float64{
			{c, s, 0, 0},
			{0, 0, c, s},
		}

		// K and M matrices
		α := o.Mdl.E * o.Mdl.A / o.L
		β := o.Mdl.GetRho() * o.Mdl.A * o.L / 6.0
		o.K = [][]float64{
			{+α * c * c, +α * c * s, -α * c * c, -α * c * s},
			{+α * c * s, +α * s * s, -α * c * s, -α * s * s},
			{-α * c * c, -α * c * s, +α * c * c, +α * c * s},
			{-α * c * s, -α * s * s, +α * c * s, +α * s * s},
		}
		o.M = [][]float64{
			{2.0 * β, 0.0, 1.0 * β, 0.0},
			{0.0, 2.0 * β, 0.0, 1.0 * β},
			{1.0 * β, 0.0, 2.0 * β, 0.0},
			{0.0, 1.0 * β, 0.0, 2.0 * β},
		}

		// return new element
		return &o
	})
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o *ElastRod) Id() int { return o.Cell.Id }

// SetEqs set equations
func (o *ElastRod) SetEqs(eqs [][]int, mixedform_eqs []int) (err error) {
	o.Umap = make([]int, o.Nu)
	for m := 0; m < 2; m++ {
		for i := 0; i < o.Ndim; i++ {
			r := i + m*o.Ndim
			o.Umap[r] = eqs[m][i]
		}
	}
	return
}

// InterpStarVars interpolates star variables to integration points
func (o *ElastRod) InterpStarVars(sol *ele.Solution) (err error) {
	// TODO: dynamics
	chk.Panic("ElastRod cannot handle dynamics yet")
	return
}

// SetEleConds set element conditions
func (o *ElastRod) SetEleConds(key string, f dbf.T, extra string) (err error) {
	if key == "g" {
		chk.Panic("ElastRod cannot handle gravity yet")
		o.Gfcn = f
	}
	return
}

// AddToRhs adds -R to global residual vector fb
func (o *ElastRod) AddToRhs(fb []float64, sol *ele.Solution) (err error) {
	for i, I := range o.Umap {
		for j, J := range o.Umap {
			fb[I] -= o.K[i][j] * sol.Y[J] // -fi
		}
	}
	return
}

// AddToKb adds element K to global Jacobian matrix Kb
func (o *ElastRod) AddToKb(Kb *la.Triplet, sol *ele.Solution, firstIt bool) (err error) {
	for i, I := range o.Umap {
		for j, J := range o.Umap {
			Kb.Put(I, J, o.K[i][j])
		}
	}
	return
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o *ElastRod) Encode(enc utl.Encoder) (err error) {
	return nil
}

// Decode decodes internal variables
func (o *ElastRod) Decode(dec utl.Decoder) (err error) {
	return nil
}

// OutIpCoords returns the coordinates of integration points
func (o *ElastRod) OutIpCoords() (C [][]float64) {
	C = utl.Alloc(1, o.Ndim) // centroid only
	for i := 0; i < o.Ndim; i++ {
		C[0][i] = (o.X[i][0] + o.X[i][1]) / 2.0 // centroid
	}
	return
}

// OutIpKeys returns the integration points' keys
func (o *ElastRod) OutIpKeys() []string {
	return []string{"sig"}
}

// OutIpVals returns the integration points' values corresponding to keys
func (o *ElastRod) OutIpVals(M *ele.IpsMap, sol *ele.Solution) {
	M.Set("sig", 0, 1, o.CalcSig(sol))
}

// specific methods /////////////////////////////////////////////////////////////////////////////////

// CalcSig computes the axial stress for given nodal displacements
func (o *ElastRod) CalcSig(sol *ele.Solution) float64 {
	for i := 0; i < 2; i++ {
		o.ua[i] = 0
		for j, J := range o.Umap {
			o.ua[i] += o.T[i][j] * sol.Y[J]
		}
	}
	εa := (o.ua[1] - o.ua[0]) / o.L // axial strain
	return o.Mdl.E * εa             // axial stress
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// Recompute re-compute matrices after dimensions or parameters are externally changed
func (o *ElastRod) Recompute(withM bool) {

	// geometry
	x0 := o.X[0][0]
	y0 := o.X[1][0]
	x1 := o.X[0][1]
	y1 := o.X[1][1]
	dx := x1 - x0
	dy := y1 - y0
	o.L = math.Sqrt(dx*dx + dy*dy)

	// global-to-local transformation matrix
	c := dx / o.L
	s := dy / o.L
	o.T[0][0] = c
	o.T[0][1] = s
	o.T[1][2] = c
	o.T[1][3] = s

	// K matrix
	α := o.Mdl.E * o.Mdl.A / o.L
	β := o.Mdl.GetRho() * o.Mdl.A * o.L / 6.0
	o.K[0][0] = +α * c * c
	o.K[0][1] = +α * c * s
	o.K[0][2] = -α * c * c
	o.K[0][3] = -α * c * s
	o.K[1][0] = +α * c * s
	o.K[1][1] = +α * s * s
	o.K[1][2] = -α * c * s
	o.K[1][3] = -α * s * s
	o.K[2][0] = -α * c * c
	o.K[2][1] = -α * c * s
	o.K[2][2] = +α * c * c
	o.K[2][3] = +α * c * s
	o.K[3][0] = -α * c * s
	o.K[3][1] = -α * s * s
	o.K[3][2] = +α * c * s
	o.K[3][3] = +α * s * s

	// M matrix
	if withM {
		o.M[0][0] = 2.0 * β
		o.M[0][2] = 1.0 * β
		o.M[1][1] = 2.0 * β
		o.M[1][3] = 1.0 * β
		o.M[2][0] = 1.0 * β
		o.M[2][2] = 2.0 * β
		o.M[3][1] = 1.0 * β
		o.M[3][3] = 2.0 * β
	}
}
