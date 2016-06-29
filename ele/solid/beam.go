// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package solid

import (
	"math"

	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mdl/sld"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

// Beam represents a structural beam element (Euler-Bernoulli, linear elastic)
//
//  2D    y1     y2 is out-of-plane
//         ^
//         | qnL          qn            qnR     Props:    Nodes:
//         o-------------------------------o     E, A      0 and 1
//         |                               |     I22
//         |                               |
//       (y2)-----------------------------(1)------> y0
//
//  3D                    ,o--------o    ,y0
//                      ,' |     ,' |  ,'
//        y1          ,'       ,'   |,'
//         ^        ,'q1     ,'    ,|
//         |      ,'  V    ,'    ,  |
//         |    ,'       ,'    ,    |
//         |  ,'       ,'  | ,      |
//         |,'       ,'   (1) - - - o   -   -  (2)
//         o--------o    ,        ,'
//         |        |  ,  <q2   ,'    Props:          Nodes:
//         |        |,        ,'       E, G, A         0, 1, 2
//         |       ,|       ,'         I22 ~ Imax      where node (2) is a point located on plane
//         |     ,  |     ,'           I11 ~ Imin      y0-y2 and non-colinear to (0) and (1).
//         |   ,    |   ,'             Jtt             Node (2) doest not have any DOF
//         | ,      | ,'
//        (0)-------o' --------> y2
//
type Beam struct {

	// basic data
	Cell *inp.Cell   // the cell structure
	X    [][]float64 // matrix of nodal coordinates [ndim][nnode]
	P02  []float64   // [3] point defining y0-y2 plane (from X matrix or computed here for horizontal/vertical beams)
	Nu   int         // total number of unknowns
	Ndim int         // space dimension

	// parameters and properties
	Mdl *sld.OnedLinElast // material model with: E, G, A, I22, I11, Jtt and Rho
	L   float64           // (derived) length of beam

	// for output
	Nstations int // number of points along beam to generate bending moment / shear force diagrams

	// variables for dynamics
	Gfcn fun.Func // gravity function

	// unit vectors aligned with beam element
	e0 []float64 // [3] unit vector aligned with y0-axis
	e1 []float64 // [3] unit vector aligned with y1-axis
	e2 []float64 // [3] unit vector aligned with y2-axis

	// vectors and matrices
	T   [][]float64 // global-to-local transformation matrix [nnode*ndim][nnode*ndim]
	Kl  [][]float64 // local K matrix
	K   [][]float64 // global K matrix
	Ml  [][]float64 // local M matrices
	M   [][]float64 // global M matrices
	Rus []float64   // residual: Rus = fi - fx

	// problem variables
	Umap []int    // assembly map (location array/element equations)
	Hasq bool     // has distributed loads
	QnL  fun.Func // distributed normal load functions: left
	QnR  fun.Func // distributed normal load functions: right
	Qt   fun.Func // distributed tangential load
	Q1   fun.Func // 3D: load on plane s-t
	Q2   fun.Func // 3D: load on plane r-t

	// scratchpad. computed @ each ip
	grav []float64 // [ndim] gravity vector
	fi   []float64 // [nu] internal forces
	ue   []float64 // local u vector
	ua   []float64 // [6] u aligned with beam system
	ζe   []float64 // local ζ* vector
	fxl  []float64 // local external force vector
}

// register element
func init() {

	// information allocator
	ele.SetInfoFunc("beam", func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *ele.Info {

		// new info
		var info ele.Info

		// solution variables
		ykeys := []string{"ux", "uy", "rz"}
		if sim.Ndim == 3 {
			ykeys = []string{"ux", "uy", "uz", "rx", "ry", "rz"}
		}
		nverts := len(cell.Verts)
		info.Dofs = make([][]string, nverts)
		for m := 0; m < 2; m++ {
			info.Dofs[m] = ykeys
		}

		// maps
		info.Y2F = map[string]string{"ux": "fx", "uy": "fy", "uz": "fz", "rx": "mx", "ry": "my", "rz": "mz"}

		// t1 and t2 variables
		info.T2vars = ykeys
		return &info
	})

	// element allocator
	ele.SetAllocator("beam", func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) ele.Element {

		// basic data
		var o Beam
		o.Cell = cell
		o.X = x
		o.P02 = []float64{0, 0, 1}
		o.Ndim = len(x)
		ndof := 3 * (o.Ndim - 1)
		o.Nu = 2 * ndof

		// model
		mat := sim.MatModels.Get(edat.Mat)
		if mat == nil {
			chk.Panic("cannot find material %q for beam {tag=%d, id=%d}\n", edat.Mat, cell.Tag, cell.Id)
		}
		o.Mdl = mat.Sld.(*sld.OnedLinElast)

		// check
		ϵp := 1e-9
		if o.Mdl.E < ϵp || o.Mdl.A < ϵp || o.Mdl.I22 < ϵp || o.Mdl.GetRho() < ϵp {
			chk.Panic("E, A, I22 and rho parameters must be all positive")
		}
		if o.Ndim == 3 {
			if o.Mdl.G < ϵp || o.Mdl.I11 < ϵp || o.Mdl.Jtt < ϵp {
				chk.Panic("G, I11, Jtt parameters must be all positive")
			}
		}

		// for output
		o.Nstations = 11
		if s_nsta, found := io.Keycode(edat.Extra, "nsta"); found {
			o.Nstations = io.Atoi(s_nsta)
		}

		// unit vectors aligned with beam element
		o.e0 = make([]float64, 3)
		o.e1 = make([]float64, 3)
		o.e2 = make([]float64, 3)

		// vectors and matrices
		o.T = la.MatAlloc(o.Nu, o.Nu)
		o.Kl = la.MatAlloc(o.Nu, o.Nu)
		o.K = la.MatAlloc(o.Nu, o.Nu)
		if !sim.Data.Steady {
			o.Ml = la.MatAlloc(o.Nu, o.Nu)
			o.M = la.MatAlloc(o.Nu, o.Nu)
		}
		o.ue = make([]float64, o.Nu)
		o.ua = make([]float64, o.Nu)
		o.ζe = make([]float64, o.Nu)
		o.fxl = make([]float64, o.Nu)
		o.Rus = make([]float64, o.Nu)

		// compute K and M
		o.Recompute(!sim.Data.Steady)

		// scratchpad. computed @ each ip
		o.grav = make([]float64, o.Ndim)
		o.fi = make([]float64, o.Nu)

		// return new element
		return &o
	})
}

// Id returns the cell Id
func (o *Beam) Id() int { return o.Cell.Id }

// SetEqs set equations [2][?]. Format of eqs == format of info.Dofs
func (o *Beam) SetEqs(eqs [][]int, mixedform_eqs []int) (err error) {
	ndof := 3 * (o.Ndim - 1)
	o.Umap = make([]int, o.Nu)
	for m := 0; m < 2; m++ {
		for i := 0; i < ndof; i++ {
			r := i + m*ndof
			o.Umap[r] = eqs[m][i]
		}
	}
	return
}

// SetEleConds set element conditions
func (o *Beam) SetEleConds(key string, f fun.Func, extra string) (err error) {

	// gravity
	if key == "g" {
		o.Gfcn = f
		return
	}

	// distributed loads
	switch key {
	case "qn":
		o.Hasq, o.QnL, o.QnR = true, f, f
	case "qnL":
		o.Hasq, o.QnL = true, f
	case "qnR":
		o.Hasq, o.QnR = true, f
	case "qt":
		o.Hasq, o.Qt = true, f
	case "q1":
		o.Hasq, o.Q1 = true, f
	case "q2":
		o.Hasq, o.Q2 = true, f
	default:
		return chk.Err("cannot handle boundary condition named %q", key)
	}
	return
}

// InterpStarVars interpolates star variables to integration points
func (o *Beam) InterpStarVars(sol *ele.Solution) (err error) {
	for i, I := range o.Umap {
		o.ζe[i] = sol.Zet[I]
	}
	return
}

// adds -R to global residual vector fb
func (o *Beam) AddToRhs(fb []float64, sol *ele.Solution) (err error) {

	// node displacements
	for i, I := range o.Umap {
		o.ue[i] = sol.Y[I]
	}

	// steady/dynamics
	if sol.Steady {
		la.MatVecMul(o.fi, 1, o.K, o.ue)
	} else {
		α1 := sol.DynCfs.GetAlp1()
		for i := 0; i < o.Nu; i++ {
			o.fi[i] = 0
			for j := 0; j < o.Nu; j++ {
				o.fi[i] += o.M[i][j]*(α1*o.ue[j]-o.ζe[j]) + o.K[i][j]*o.ue[j]
			}
		}
	}

	// distributed loads
	if o.Hasq {
		l := o.L
		ll := l * l
		qnL, qnR, qt, q1, q2 := o.calc_loads(sol.T)
		if o.Ndim == 2 {
			o.fxl[0] = qt * l / 2.0
			o.fxl[1] = l * (7.0*qnL + 3.0*qnR) / 20.0
			o.fxl[2] = ll * (3.0*qnL + 2.0*qnR) / 60.0
			o.fxl[3] = qt * l / 2.0
			o.fxl[4] = l * (3.0*qnL + 7.0*qnR) / 20.0
			o.fxl[5] = -ll * (2.0*qnL + 3.0*qnR) / 60.0
		} else {
			o.fxl[1] = l * q1 / 2.0
			o.fxl[2] = l * q2 / 2.0
			o.fxl[4] = -ll * q2 / 12.0
			o.fxl[5] = ll * q1 / 12.0
			o.fxl[7] = l * q1 / 2.0
			o.fxl[8] = l * q2 / 2.0
			o.fxl[10] = ll * q2 / 12.0
			o.fxl[11] = -ll * q1 / 12.0
		}
		la.MatTrVecMulAdd(o.fi, -1.0, o.T, o.fxl) // Rus -= fx; fx = trans(T) * fxl
	}

	// add to fb
	for i, I := range o.Umap {
		fb[I] -= o.fi[i]
	}
	return
}

// adds element K to global Jacobian matrix Kb
func (o *Beam) AddToKb(Kb *la.Triplet, sol *ele.Solution, firstIt bool) (err error) {
	if sol.Steady {
		for i, I := range o.Umap {
			for j, J := range o.Umap {
				Kb.Put(I, J, o.K[i][j])
			}
		}
		return
	}
	α1 := sol.DynCfs.GetAlp1()
	for i, I := range o.Umap {
		for j, J := range o.Umap {
			Kb.Put(I, J, o.M[i][j]*α1+o.K[i][j])
		}
	}
	return
}

// Encode encodes internal variables
func (o *Beam) Encode(enc utl.Encoder) (err error) {
	return
}

// Decode decodes internal variables
func (o *Beam) Decode(dec utl.Decoder) (err error) {
	return
}

// OutIpCoords returns the coordinates of integration points
func (o *Beam) OutIpCoords() (C [][]float64) {
	C = make([][]float64, o.Nstations)
	dξ := 1.0 / float64(o.Nstations-1)
	for i := 0; i < o.Nstations; i++ {
		ξ := float64(i) * dξ
		C[i] = make([]float64, o.Ndim)
		for j := 0; j < o.Ndim; j++ {
			C[i][j] = (1.0-ξ)*o.X[j][0] + ξ*o.X[j][1]
		}
	}
	return
}

// OutIpKeys returns the integration points' keys
func (o *Beam) OutIpKeys() []string {
	if o.Ndim == 3 {
		return []string{"M22", "M11", "T00"}
	}
	return []string{"M22"}
}

// OutIpVals returns the integration points' values corresponding to keys
func (o *Beam) OutIpVals(M *ele.IpsMap, sol *ele.Solution) {
	unused := 0
	dξ := 1.0 / float64(o.Nstations-1)
	for i := 0; i < o.Nstations; i++ {
		ξ := float64(i) * dξ
		if o.Ndim == 3 {
			M22, M11, T00 := o.CalcMoment3d(sol, ξ, unused)
			M.Set("M22", i, o.Nstations, M22[0])
			M.Set("M11", i, o.Nstations, M11[0])
			M.Set("T00", i, o.Nstations, T00[0])
		} else {
			M22 := o.CalcMoment2d(sol, ξ, unused)
			M.Set("M22", i, o.Nstations, M22[0])
		}
	}
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// Recompute re-compute matrices after dimensions or parameters are externally changed
func (o *Beam) Recompute(withM bool) {

	// 3D
	if o.Ndim == 3 {

		// point defining y0-y2 plane
		if len(o.X[0]) == 3 { // point given
			for i := 0; i < o.Ndim; i++ {
				o.P02[i] = o.X[i][2]
			}
		} else {
			dx := make([]float64, 3)
			for i := 0; i < 3; i++ {
				dx[i] = o.X[i][1] - o.X[i][0]
			}
			tol := 1e-5 // tolerance to find horizontal/vertical beams
			switch {

			// vertical (parallel to z)
			case math.Abs(dx[0]) < tol && math.Abs(dx[1]) < tol:
				δ := 0.1 * dx[2] // + if 0->1 is going up
				o.P02[0], o.P02[1], o.P02[2] = o.X[0][0]+δ, o.X[1][0], o.X[2][0]

			// horizontal (perpendicular to z)
			case math.Abs(dx[2]) < tol:
				o.e0[0], o.e0[1], o.e0[2] = dx[0], dx[1], 0
				o.e1[0], o.e1[1], o.e1[2] = 0, 0, 1
				utl.Cross3d(o.e2, o.e0, o.e1) // e2 := e0 cross e1
				l0 := la.VecNorm(o.e0)
				l2 := la.VecNorm(o.e2)
				δ := 0.1 * l0 / l2
				o.P02[0], o.P02[1], o.P02[2] = o.X[0][0]+δ*o.e2[0], o.X[1][0]+δ*o.e2[1], o.X[2][0]

			default:
				chk.Panic("Beam: can only compute P02 vertex for vertical and horizontal beams")
			}
		}

		// auxiliary vector
		o.L = 0.0
		v02 := make([]float64, o.Ndim)
		for i := 0; i < o.Ndim; i++ {
			o.e0[i] = o.X[i][1] - o.X[i][0]
			v02[i] = o.P02[i] - o.X[i][0]
			o.L += o.e0[i] * o.e0[i]
		}
		o.L = math.Sqrt(o.L)
		utl.Cross3d(o.e1, v02, o.e0) // e1 := v02 cross e0

		// unit vectors aligned with beam element
		nrm1 := la.VecNorm(o.e1)
		for i := 0; i < o.Ndim; i++ {
			o.e0[i] = o.e0[i] / o.L
			o.e1[i] = o.e1[i] / nrm1
		}
		utl.Cross3d(o.e2, o.e0, o.e1) // e2 := e0 cross e1

		// global to local transformation matrix
		for k := 0; k < 4; k++ {
			o.T[3*k+0][3*k+0], o.T[3*k+0][3*k+1], o.T[3*k+0][3*k+2] = o.e0[0], o.e0[1], o.e0[2]
			o.T[3*k+1][3*k+0], o.T[3*k+1][3*k+1], o.T[3*k+1][3*k+2] = o.e1[0], o.e1[1], o.e1[2]
			o.T[3*k+2][3*k+0], o.T[3*k+2][3*k+1], o.T[3*k+2][3*k+2] = o.e2[0], o.e2[1], o.e2[2]
		}

		// constants
		EIr := o.Mdl.E * o.Mdl.I22
		EIs := o.Mdl.E * o.Mdl.I11
		GJ := o.Mdl.G * o.Mdl.Jtt
		EA := o.Mdl.E * o.Mdl.A
		l := o.L
		ll := l * l
		lll := l * ll

		// stiffness matrix in local system
		o.Kl[0][0] = EA / l  //  0
		o.Kl[0][6] = -EA / l //  1

		o.Kl[1][1] = 12.0 * EIr / lll  //  2
		o.Kl[1][5] = 6.0 * EIr / ll    //  3
		o.Kl[1][7] = -12.0 * EIr / lll //  4
		o.Kl[1][11] = 6.0 * EIr / ll   //  5

		o.Kl[2][2] = 12.0 * EIs / lll  //  6
		o.Kl[2][4] = -6.0 * EIs / ll   //  7
		o.Kl[2][8] = -12.0 * EIs / lll //  8
		o.Kl[2][10] = -6.0 * EIs / ll  //  9

		o.Kl[3][3] = GJ / l  // 10
		o.Kl[3][9] = -GJ / l // 11

		o.Kl[4][2] = -6.0 * EIs / ll // 12
		o.Kl[4][4] = 4.0 * EIs / l   // 13
		o.Kl[4][8] = 6.0 * EIs / ll  // 14
		o.Kl[4][10] = 2.0 * EIs / l  // 15

		o.Kl[5][1] = 6.0 * EIr / ll  // 16
		o.Kl[5][5] = 4.0 * EIr / l   // 17
		o.Kl[5][7] = -6.0 * EIr / ll // 18
		o.Kl[5][11] = 2.0 * EIr / l  // 19

		o.Kl[6][0] = -EA / l // 20
		o.Kl[6][6] = EA / l  // 21

		o.Kl[7][1] = -12.0 * EIr / lll // 22
		o.Kl[7][5] = -6.0 * EIr / ll   // 23
		o.Kl[7][7] = 12.0 * EIr / lll  // 24
		o.Kl[7][11] = -6.0 * EIr / ll  // 25

		o.Kl[8][2] = -12.0 * EIs / lll // 26
		o.Kl[8][4] = 6.0 * EIs / ll    // 27
		o.Kl[8][8] = 12.0 * EIs / lll  // 28
		o.Kl[8][10] = 6.0 * EIs / ll   // 29

		o.Kl[9][3] = -GJ / l // 30
		o.Kl[9][9] = GJ / l  // 31

		o.Kl[10][2] = -6.0 * EIs / ll // 32
		o.Kl[10][4] = 2.0 * EIs / l   // 33
		o.Kl[10][8] = 6.0 * EIs / ll  // 34
		o.Kl[10][10] = 4.0 * EIs / l  // 35

		o.Kl[11][1] = 6.0 * EIr / ll  // 36
		o.Kl[11][5] = 2.0 * EIr / l   // 37
		o.Kl[11][7] = -6.0 * EIr / ll // 38
		o.Kl[11][11] = 4.0 * EIr / l  // 39

		// stiffness matrix in global system
		la.MatTrMul3(o.K, 1, o.T, o.Kl, o.T) // K := 1 * trans(T) * Kl * T

		// mass matrix
		if withM {
			chk.Panic("mass matrix is not available for 3D beams yet")
		}
		return
	}

	// T
	dx := o.X[0][1] - o.X[0][0]
	dy := o.X[1][1] - o.X[1][0]
	l := math.Sqrt(dx*dx + dy*dy)
	o.L = l
	c := dx / l
	s := dy / l
	o.T[0][0] = c
	o.T[0][1] = s
	o.T[1][0] = -s
	o.T[1][1] = c
	o.T[2][2] = 1
	o.T[3][3] = c
	o.T[3][4] = s
	o.T[4][3] = -s
	o.T[4][4] = c
	o.T[5][5] = 1

	// unit vectors aligned with beam element
	o.e0[0], o.e0[1] = c, s
	o.e1[0], o.e1[1] = -s, c
	o.e2[2] = 1

	// aux vars
	ll := l * l
	m := o.Mdl.E * o.Mdl.A / l
	n := o.Mdl.E * o.Mdl.I22 / (ll * l)

	// K
	o.Kl[0][0] = m
	o.Kl[0][3] = -m
	o.Kl[1][1] = 12 * n
	o.Kl[1][2] = 6 * l * n
	o.Kl[1][4] = -12 * n
	o.Kl[1][5] = 6 * l * n
	o.Kl[2][1] = 6 * l * n
	o.Kl[2][2] = 4 * ll * n
	o.Kl[2][4] = -6 * l * n
	o.Kl[2][5] = 2 * ll * n
	o.Kl[3][0] = -m
	o.Kl[3][3] = m
	o.Kl[4][1] = -12 * n
	o.Kl[4][2] = -6 * l * n
	o.Kl[4][4] = 12 * n
	o.Kl[4][5] = -6 * l * n
	o.Kl[5][1] = 6 * l * n
	o.Kl[5][2] = 2 * ll * n
	o.Kl[5][4] = -6 * l * n
	o.Kl[5][5] = 4 * ll * n
	la.MatTrMul3(o.K, 1, o.T, o.Kl, o.T) // K := 1 * trans(T) * Kl * T

	// M
	if withM {
		m = o.Mdl.GetRho() * o.Mdl.A * l / 420.0
		o.Ml[0][0] = 140.0 * m
		o.Ml[0][3] = 70.0 * m
		o.Ml[1][1] = 156.0 * m
		o.Ml[1][2] = 22.0 * l * m
		o.Ml[1][4] = 54.0 * m
		o.Ml[1][5] = -13.0 * l * m
		o.Ml[2][1] = 22.0 * l * m
		o.Ml[2][2] = 4.0 * ll * m
		o.Ml[2][4] = 13.0 * l * m
		o.Ml[2][5] = -3.0 * ll * m
		o.Ml[3][0] = 70.0 * m
		o.Ml[3][3] = 140.0 * m
		o.Ml[4][1] = 54.0 * m
		o.Ml[4][2] = 13.0 * l * m
		o.Ml[4][4] = 156.0 * m
		o.Ml[4][5] = -22.0 * l * m
		o.Ml[5][1] = -13.0 * l * m
		o.Ml[5][2] = -3.0 * ll * m
		o.Ml[5][4] = -22.0 * l * m
		o.Ml[5][5] = 4.0 * ll * m
		la.MatTrMul3(o.M, 1, o.T, o.Ml, o.T) // M := 1 * trans(T) * Ml * T
	}
}

// calc_loads computes applied distributed loads at given time
func (o *Beam) calc_loads(time float64) (qnL, qnR, qt, q1, q2 float64) {
	if o.QnL != nil {
		qnL = o.QnL.F(time, nil)
	}
	if o.QnR != nil {
		qnR = o.QnR.F(time, nil)
	}
	if o.Qt != nil {
		qt = o.Qt.F(time, nil)
	}
	if o.Q1 != nil {
		q1 = o.Q1.F(time, nil)
	}
	if o.Q2 != nil {
		q2 = o.Q2.F(time, nil)
	}
	return
}

// bending moments and forces ///////////////////////////////////////////////////////////////////////

// calc_ua computes local (aligned) displacements
func (o *Beam) calc_ua(sol *ele.Solution) {
	for i := 0; i < o.Nu; i++ {
		o.ua[i] = 0
		for j, J := range o.Umap {
			o.ua[i] += o.T[i][j] * sol.Y[J]
		}
	}
}

// CalcMoment3d calculates moments along 3D beam
//  Input:
//   ξ         -- natural coordinate along bar   0 ≤ ξ ≤ 1
//   nstations -- compute many values; otherwise, if nstations<2, compute @ s
//  Output:
//   M22 -- bending moment about y2-axis
//   M11 -- bending moment about y1-axis
//   T00 -- twisting moment around y0-axis
func (o *Beam) CalcMoment3d(sol *ele.Solution, ξ float64, nstations int) (M22, M11, T00 []float64) {
	o.calc_ua(sol)
	if nstations < 2 {
		mrr, mss, mtt := o.calc_moment3d_after_ua(sol.T, ξ)
		M22, M11, T00 = []float64{mrr}, []float64{mss}, []float64{mtt}
		return
	}
	M22 = make([]float64, nstations)
	M11 = make([]float64, nstations)
	T00 = make([]float64, nstations)
	dξ := 1.0 / float64(nstations-1)
	for i := 0; i < nstations; i++ {
		M22[i], M11[i], T00[i] = o.calc_moment3d_after_ua(sol.T, float64(i)*dξ)
	}
	return
}

// CalcMoment2d calculates bending moment along 2D beam
//  Input:
//   ξ         -- natural coordinate along bar   0 ≤ ξ ≤ 1
//   nstations -- compute many values; otherwise, if nstations<2, compute @ s
//  Output:
//   M22 -- bending moment @ stations or s
func (o *Beam) CalcMoment2d(sol *ele.Solution, ξ float64, nstations int) (M22 []float64) {
	o.calc_ua(sol)
	if nstations < 2 {
		m := o.calc_moment2d_after_ua(sol.T, ξ)
		M22 = []float64{m}
		return
	}
	M22 = make([]float64, nstations)
	dξ := 1.0 / float64(nstations-1)
	for i := 0; i < nstations; i++ {
		M22[i] = o.calc_moment2d_after_ua(sol.T, float64(i)*dξ)
	}
	return
}

// CalcShearForce2d calculates shear force for 2D beam
//  Input:
//   ξ         -- natural coordinate along bar   0 ≤ ξ ≤ 1
//   nstations -- compute many values; otherwise, if nstations<2, compute @ s
//  Output:
//   V1 -- shear force @ stations or s
func (o *Beam) CalcShearForce2d(sol *ele.Solution, ξ float64, nstations int) (V1 []float64) {
	o.calc_ua(sol)
	if nstations < 2 {
		v := o.calc_shearforce2d_after_ua(sol.T, ξ)
		V1 = []float64{v}
		return
	}
	V1 = make([]float64, nstations)
	dξ := 1.0 / float64(nstations-1)
	for i := 0; i < nstations; i++ {
		V1[i] = o.calc_shearforce2d_after_ua(sol.T, float64(i)*dξ)
	}
	return
}

// calc_bendingmom3d_after_ua calculates bending moments and torque (3D) @ station ξ in [0, 1]
func (o *Beam) calc_moment3d_after_ua(time, ξ float64) (M22, M11, T00 float64) {

	// auxiliary variables
	τ := ξ * o.L
	τ2 := τ * τ
	l := o.L
	ll := l * l
	lll := ll * l

	// constants and loads
	EIr := o.Mdl.E * o.Mdl.I22
	EIs := o.Mdl.E * o.Mdl.I11
	GJ := o.Mdl.G * o.Mdl.Jtt
	_, _, _, q1, q2 := o.calc_loads(time)

	// displacements and rotations
	//u := []float64{o.ua[0], o.ua[6]}
	v := []float64{o.ua[1], o.ua[5], o.ua[7], o.ua[11]}
	w := []float64{o.ua[2], o.ua[4], o.ua[8], o.ua[10]}
	θ := []float64{o.ua[3], o.ua[9]}

	// second derivatives of shape functions
	dnv2 := []float64{12.0*τ/lll - 6.0/ll, 6.0*τ/ll - 4.0/l, 6.0/ll - 12.0*τ/lll, 6.0*τ/ll - 2.0/l}
	dnw2 := []float64{dnv2[0], -dnv2[1], dnv2[2], -dnv2[3]}

	// bending moments
	M22 = +q1 * (ll - 6*τ*l + 6*τ2) / 12.0 //   EIr*dnv2*v
	M11 = -q2 * (ll - 6*τ*l + 6*τ2) / 12.0 //  -EIs*dnw2*w
	for i := 0; i < len(v); i++ {
		M22 += EIr * dnv2[i] * v[i]
		M11 -= EIs * dnw2[i] * w[i]
	}
	T00 = GJ * (θ[1] - θ[0]) / l
	return
}

// calc_bendingmom2d_after_ua calculates bending moment (2D) @ station ξ in [0, 1]
func (o *Beam) calc_moment2d_after_ua(time, ξ float64) (M22 float64) {

	// auxiliary variables
	τ := ξ * o.L
	l := o.L
	ll := l * l
	lll := ll * l

	// bending moment
	M22 = o.Mdl.E * o.Mdl.I22 * (o.ua[1]*((12.0*τ)/lll-6.0/ll) + o.ua[2]*((6.0*τ)/ll-4.0/l) + o.ua[4]*(6.0/ll-(12.0*τ)/lll) + o.ua[5]*((6.0*τ)/ll-2.0/l))

	// corrections due to applied loads
	if o.Hasq {
		qnL, qnR, _, _, _ := o.calc_loads(time)
		τ2 := τ * τ
		τ3 := τ2 * τ
		M22 += (2.0*qnR*lll + 3.0*qnL*lll - 9.0*qnR*τ*ll - 21.0*qnL*τ*ll + 30.0*qnL*τ2*l + 10.0*qnR*τ3 - 10.0*qnL*τ3) / (60.0 * l)
		if qnL > 0.0 {
			M22 = -M22 // swap the sign of M
		}
	}
	return
}

// calc_shearforce2d_after_ua calculates shear force (2D) @ station ξ in [0, 1]
func (o *Beam) calc_shearforce2d_after_ua(time, ξ float64) (V1 float64) {

	// auxiliary variables
	τ := ξ * o.L
	l := o.L
	ll := l * l
	lll := ll * l

	// shear force
	V1 = o.Mdl.E * o.Mdl.I22 * ((12.0*o.ua[1])/lll + (6.0*o.ua[2])/ll - (12.0*o.ua[4])/lll + (6.0*o.ua[5])/ll)

	// corrections due to applied loads
	if o.Hasq {
		qnL, qnR, _, _, _ := o.calc_loads(time)
		τ2 := τ * τ
		V1 += -(3.0*qnR*ll + 7.0*qnL*ll - 20.0*qnL*τ*l - 10.0*qnR*τ2 + 10.0*qnL*τ2) / (20.0 * l)
	}
	return
}

// plot diagrams ////////////////////////////////////////////////////////////////////////////////////

// PlotDiagMoment plots bending moment diagram
//  Input:
//   M        -- moment along stations
//   withtext -- show bending moment values
//   numfmt   -- number format for values. use "" to chose default one
//   tolM     -- tolerance to clip absolute values of M
//   sf       -- scaling factor
func (o *Beam) PlotDiagMoment(M []float64, withtext bool, numfmt string, tolM, sf float64) {

	// number of stations
	nstations := len(M)
	ds := 1.0 / float64(nstations-1)

	// nodes
	var xa, xb []float64
	var u []float64 // out-of-pane vector
	if o.Ndim == 2 {
		xa = []float64{o.X[0][0], o.X[1][0], 0}
		xb = []float64{o.X[0][1], o.X[1][1], 0}
		u = []float64{0, 0, 1}
	} else {
		chk.Panic("TODO: 3D beam diagram")
	}

	// unit vector along beam
	v := make([]float64, 3)
	sum := 0.0
	for j := 0; j < o.Ndim; j++ {
		v[j] = xb[j] - xa[j]
		sum += v[j] * v[j]
	}
	sum = math.Sqrt(sum)
	for j := 0; j < o.Ndim; j++ {
		v[j] /= sum
	}

	// unit normal
	n := make([]float64, 3) // normal
	utl.Cross3d(n, u, v)    // n := u cross v

	// auxiliary vectors
	x := make([]float64, o.Ndim) // station
	m := make([]float64, o.Ndim) // vector pointing to other side
	c := make([]float64, o.Ndim) // centre
	imin, imax := utl.DblArgMinMax(M)

	// draw text function
	draw_text := func(mom float64) {
		if math.Abs(mom) > tolM {
			α := math.Atan2(-n[1], -n[0]) * 180.0 / math.Pi
			str := io.Sf("%g", mom)
			if numfmt != "" {
				str = io.Sf(numfmt, mom)
			} else {
				if len(str) > 10 {
					str = io.Sf("%.8f", mom) // truncate number
					str = io.Sf("%g", io.Atof(str))
				}
			}
			plt.Text(c[0], c[1], str, io.Sf("ha='center', size=7, rotation=%g, clip_on=0", α))
		}
	}

	// draw
	pts := utl.DblsAlloc(nstations, 2)
	xx, yy := make([]float64, 2), make([]float64, 2)
	for i := 0; i < nstations; i++ {

		// station
		s := float64(i) * ds
		for j := 0; j < o.Ndim; j++ {
			x[j] = (1.0-s)*o.X[j][0] + s*o.X[j][1]
		}

		// auxiliary vectors
		for j := 0; j < o.Ndim; j++ {
			m[j] = x[j] - sf*M[i]*n[j]
			c[j] = (x[j] + m[j]) / 2.0
		}

		// points on diagram
		pts[i][0], pts[i][1] = m[0], m[1]
		xx[0], xx[1] = x[0], m[0]
		yy[0], yy[1] = x[1], m[1]

		// draw
		clr, lw := "#919191", 1.0
		if i == imin || i == imax {
			lw = 2
			if M[i] < 0 {
				clr = "#9f0000"
			} else {
				clr = "#109f24"
			}
		}
		plt.Plot(xx, yy, io.Sf("'-', color='%s', lw=%g, clip_on=0", clr, lw))
		if withtext {
			if i == imin || i == imax { // draw text @ min/max
				draw_text(M[i])
			} else {
				if i == 0 || i == nstations-1 { // draw text @ extremities
					draw_text(M[i])
				}
			}
		}
	}

	// draw polyline
	plt.DrawPolyline(pts, &plt.Sty{Ec: "k", Fc: "none", Lw: 1}, "")
}
