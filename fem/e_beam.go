// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"

	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

// Beam represents a structural beam element (Euler-Bernoulli, linear elastic)
//
//  2D     s     r is out-of-plane              Note: r,s,t not under any right-hand-type rule
//         ^
//         |                                    Props:  Nodes:
//         o-------------------------------o     E, A    0 and 1
//         |                               |     Irr
//         |                               |
//        (0)-----------------------------(1)------> t
//
//  3D                    ,o--------o    ,t
//                      ,' |     ,' |  ,'
//         s          ,'       ,'   |,'
//         ^        ,'       ,'    ,|
//         |      ,'       ,'    ,  |
//         |    ,'       ,'    ,    |
//         |  ,'       ,'  | ,      |
//         |,'       ,'   (1) - - - o   -   -  (2)
//         o--------o    ,        ,'
//         |        |  ,        ,'    Props:       Nodes:
//         |        |,        ,'       E, G, A      0, 1, 2
//         |       ,|       ,'         Irr = Imax   where node (2) is a point located on plane r-t
//         |     ,  |     ,'           Iss = Imin   and non-colinear to (0) and (1)
//         |   ,    |   ,'             Jtt          Node (2) doest not have any DOF
//         | ,      | ,'
//        (0)-------o' --------> r
//
type Beam struct {

	// basic data
	Cell *inp.Cell   // the cell structure
	X    [][]float64 // matrix of nodal coordinates [ndim][nnode]
	Nu   int         // total number of unknowns
	Ndim int         // space dimension

	// parameters and properties
	E   float64 // Young's modulus
	G   float64 // shear modulus
	A   float64 // cross-sectional area
	Irr float64 // moment of inertia of cross section about r-axis (maximum principal inertia)
	Iss float64 // moment of inertia of cross section about s-axis (minimum principal inertia)
	Jtt float64 // torsional constant
	L   float64 // (derived) length of beam

	// for output
	Nstations int // number of points along beam to generate bending moment / shear force diagrams

	// variables for dynamics
	Rho  float64  // density of solids
	Gfcn fun.Func // gravity function

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
	Qs   fun.Func // 3D: load on plane s-t
	Qr   fun.Func // 3D: load on plane r-t

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
	infogetters["beam"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *Info {

		// new info
		var info Info

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
	}

	// element allocator
	eallocators["beam"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) Elem {

		// basic data
		var o Beam
		o.Cell = cell
		o.X = x
		o.Ndim = len(x)
		ndof := 3 * (o.Ndim - 1)
		o.Nu = 2 * ndof

		// parameters
		matdata := sim.MatParams.Get(edat.Mat)
		if matdata == nil {
			return nil
		}
		for _, p := range matdata.Prms {
			switch p.N {
			case "E":
				o.E = p.V
			case "G":
				o.G = p.V
			case "A":
				o.A = p.V
			case "Irr", "Imax":
				o.Irr = p.V
			case "Iss", "Imin":
				o.Iss = p.V
			case "Jtt":
				o.Jtt = p.V
			case "rho":
				o.Rho = p.V
			}
		}
		ϵp := 1e-9
		if o.E < ϵp || o.A < ϵp || o.Irr < ϵp || o.Rho < ϵp {
			chk.Panic("E, A, Irr and rho parameters must be all positive")
		}
		if o.Ndim == 3 {
			if o.G < ϵp || o.Iss < ϵp || o.Jtt < ϵp {
				chk.Panic("G, Iss, Jtt parameters must be all positive")
			}
		}

		// for output
		o.Nstations = 11
		if s_nsta, found := io.Keycode(edat.Extra, "nsta"); found {
			o.Nstations = io.Atoi(s_nsta)
		}

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
	}
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
	case "qs":
		o.Hasq, o.Qs = true, f
	case "qr":
		o.Hasq, o.Qr = true, f
	default:
		return chk.Err("cannot handle boundary condition named %q", key)
	}
	return
}

// InterpStarVars interpolates star variables to integration points
func (o *Beam) InterpStarVars(sol *Solution) (err error) {
	for i, I := range o.Umap {
		o.ζe[i] = sol.Zet[I]
	}
	return
}

// adds -R to global residual vector fb
func (o *Beam) AddToRhs(fb []float64, sol *Solution) (err error) {

	// node displacements
	for i, I := range o.Umap {
		o.ue[i] = sol.Y[I]
	}

	// steady/dynamics
	if sol.Steady {
		la.MatVecMul(o.fi, 1, o.K, o.ue)
	} else {
		α1 := sol.DynCfs.α1
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
		qnL, qnR, qt, qs, qr := o.calc_loads(sol.T)
		if o.Ndim == 2 {
			o.fxl[0] = qt * l / 2.0
			o.fxl[1] = l * (7.0*qnL + 3.0*qnR) / 20.0
			o.fxl[2] = ll * (3.0*qnL + 2.0*qnR) / 60.0
			o.fxl[3] = qt * l / 2.0
			o.fxl[4] = l * (3.0*qnL + 7.0*qnR) / 20.0
			o.fxl[5] = -ll * (2.0*qnL + 3.0*qnR) / 60.0
		} else {
			o.fxl[1] = l * qs / 2.0
			o.fxl[2] = l * qr / 2.0
			o.fxl[4] = -ll * qr / 12.0
			o.fxl[5] = ll * qs / 12.0
			o.fxl[7] = l * qs / 2.0
			o.fxl[8] = l * qr / 2.0
			o.fxl[10] = ll * qr / 12.0
			o.fxl[11] = -ll * qs / 12.0
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
func (o *Beam) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (err error) {
	if sol.Steady {
		for i, I := range o.Umap {
			for j, J := range o.Umap {
				Kb.Put(I, J, o.K[i][j])
			}
		}
		return
	}
	α1 := sol.DynCfs.α1
	for i, I := range o.Umap {
		for j, J := range o.Umap {
			Kb.Put(I, J, o.M[i][j]*α1+o.K[i][j])
		}
	}
	return
}

// Update perform (tangent) update
func (o *Beam) Update(sol *Solution) (err error) {
	return
}

// Encode encodes internal variables
func (o *Beam) Encode(enc Encoder) (err error) {
	return
}

// Decode decodes internal variables
func (o *Beam) Decode(dec Decoder) (err error) {
	return
}

// OutIpsData returns data from all integration points for output
func (o *Beam) OutIpsData() (data []*OutIpData) {
	unused := 0
	ds := 1.0 / float64(o.Nstations-1)
	for i := 0; i < o.Nstations; i++ {
		s := float64(i) * ds
		x := make([]float64, o.Ndim)
		for j := 0; j < o.Ndim; j++ {
			x[j] = (1.0-s)*o.X[j][0] + s*o.X[j][1]
		}
		calc := func(sol *Solution) (vals map[string]float64) {
			vals = make(map[string]float64)
			V, M := o.CalcVandM2d(sol, s, unused)
			vals["V"] = V[0]
			vals["M"] = M[0]
			return
		}
		data = append(data, &OutIpData{o.Id(), x, calc})
	}
	return
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// Recompute re-compute matrices after dimensions or parameters are externally changed
func (o *Beam) Recompute(withM bool) {

	// 3D
	if o.Ndim == 3 {

		// auxiliary vectors
		v01, v02 := make([]float64, o.Ndim), make([]float64, o.Ndim)
		for i := 0; i < o.Ndim; i++ {
			v01[i] = o.X[i][1] - o.X[i][0]
			v02[i] = o.X[i][2] - o.X[i][0]
		}
		vs := make([]float64, o.Ndim)
		utl.Cross3d(vs, v02, v01) // vs := v02 cross v01
		l := math.Sqrt(utl.Dot3d(v01, v01))
		o.L = l
		ls := math.Sqrt(utl.Dot3d(vs, vs))
		vt := make([]float64, o.Ndim)
		for i := 0; i < o.Ndim; i++ {
			vt[i] = v01[i] / l
			vs[i] = vs[i] / ls
		}
		vr := make([]float64, o.Ndim)
		utl.Cross3d(vr, vt, vs) // vr := vt cross vs

		// global to local transformation matrix
		for k := 0; k < 4; k++ {
			o.T[3*k+0][3*k+0], o.T[3*k+0][3*k+1], o.T[3*k+0][3*k+2] = vt[0], vt[1], vt[2]
			o.T[3*k+1][3*k+0], o.T[3*k+1][3*k+1], o.T[3*k+1][3*k+2] = vs[0], vs[1], vs[2]
			o.T[3*k+2][3*k+0], o.T[3*k+2][3*k+1], o.T[3*k+2][3*k+2] = vr[0], vr[1], vr[2]
		}

		// constants
		EIr, EIs, GJ, EA := o.E*o.Irr, o.E*o.Iss, o.G*o.Jtt, o.E*o.A
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

	// aux vars
	ll := l * l
	m := o.E * o.A / l
	n := o.E * o.Irr / (ll * l)

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
		m = o.Rho * o.A * l / 420.0
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

// CalcVandM calculate shear force and bending moment @ s
//  Input:
//   s         -- natural coordinate   0 ≤ s ≤ 1
//   nstations -- compute many values; otherwise, if nstations<2, compute @ s
//  Output:
//   V -- shear force @ stations or s
//   M -- bending moment @ stations or s
func (o *Beam) CalcVandM2d(sol *Solution, s float64, nstations int) (V, M []float64) {

	// aligned displacements
	for i := 0; i < 6; i++ {
		o.ua[i] = 0
		for j, J := range o.Umap {
			o.ua[i] += o.T[i][j] * sol.Y[J]
		}
	}

	// results
	if nstations < 2 {
		v, m := o.calc_V_and_M_after_ua2d(sol.T, s)
		V, M = []float64{v}, []float64{m}
		return
	}
	V = make([]float64, nstations)
	M = make([]float64, nstations)
	ds := 1.0 / float64(nstations-1)
	for i := 0; i < nstations; i++ {
		V[i], M[i] = o.calc_V_and_M_after_ua2d(sol.T, float64(i)*ds)
	}
	return
}

func (o *Beam) calc_V_and_M_after_ua2d(time, s float64) (V, M float64) {

	// auxiliary variables
	r := s * o.L
	l := o.L
	ll := l * l
	lll := ll * l

	// shear force
	V = o.E * o.Irr * ((12.0*o.ua[1])/lll + (6.0*o.ua[2])/ll - (12.0*o.ua[4])/lll + (6.0*o.ua[5])/ll)

	// bending moment
	M = o.E * o.Irr * (o.ua[1]*((12.0*r)/lll-6.0/ll) + o.ua[2]*((6.0*r)/ll-4.0/l) + o.ua[4]*(6.0/ll-(12.0*r)/lll) + o.ua[5]*((6.0*r)/ll-2.0/l))

	// corrections due to applied loads
	if o.Hasq {
		qnL, qnR, _, _, _ := o.calc_loads(time)
		rr := r * r
		rrr := rr * r
		V += -(3.0*qnR*ll + 7.0*qnL*ll - 20.0*qnL*r*l - 10.0*qnR*rr + 10.0*qnL*rr) / (20.0 * l)
		M += (2.0*qnR*lll + 3.0*qnL*lll - 9.0*qnR*r*ll - 21.0*qnL*r*ll + 30.0*qnL*rr*l + 10.0*qnR*rrr - 10.0*qnL*rrr) / (60.0 * l)
		if qnL > 0.0 {
			M = -M // swap the sign of M
		}
	}
	return
}

func (o *Beam) calc_loads(time float64) (qnL, qnR, qt, qs, qr float64) {
	if o.QnL != nil {
		qnL = o.QnL.F(time, nil)
	}
	if o.QnR != nil {
		qnR = o.QnR.F(time, nil)
	}
	if o.Qt != nil {
		qt = o.Qt.F(time, nil)
	}
	if o.Qs != nil {
		qs = o.Qs.F(time, nil)
	}
	if o.Qr != nil {
		qr = o.Qr.F(time, nil)
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
