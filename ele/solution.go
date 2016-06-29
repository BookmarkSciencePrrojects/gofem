// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package ele

// Solution holds the solution data @ nodes.
//
//        / u \         / u \
//        |   | => y =  |   |
//  yb =  | p |         \ p / (ny x 1)
//        |   |
//        \ λ / (nyb x 1)
//
type Solution struct {

	// current state
	T      float64   // current time
	Y      []float64 // DOFs (solution variables); e.g. y = {u, p}
	Dydt   []float64 // dy/dt
	D2ydt2 []float64 // d²y/dt²

	// auxiliary
	Dt  float64   // current time increment
	ΔY  []float64 // total increment (for nonlinear solver)
	Psi []float64 // t1 star vars; e.g. ψ* = β1.p + β2.dpdt
	Zet []float64 // t2 star vars; e.g. ζ* = α1.u + α2.v + α3.a
	Chi []float64 // t2 star vars; e.g. χ* = α4.u + α5.v + α6.a
	L   []float64 // Lagrange multipliers

	// extrapolated values
	Ext map[int][]float64 // [optional] extrapolated values. nodeId => values (e.g. σ)
	Cnt map[int]int       // [optional] counter for number of additions each node of Ext

	// problem definition and constants
	Steady  bool      // [from Sim] steady simulation
	Axisym  bool      // [from Sim] axisymmetric
	Pstress bool      // [from Sim] plane-stress
	DynCfs  *DynCoefs // [from FEM] coefficients for dynamics/transient simulations
}

// Reset clear values
func (o *Solution) Reset(steady bool) {
	o.T = 0
	for i := 0; i < len(o.Y); i++ {
		o.Y[i] = 0
		o.ΔY[i] = 0
	}
	if !steady {
		for i := 0; i < len(o.Y); i++ {
			o.Psi[i] = 0
			o.Zet[i] = 0
			o.Chi[i] = 0
			o.Dydt[i] = 0
			o.D2ydt2[i] = 0
		}
	}
	for i := 0; i < len(o.L); i++ {
		o.L[i] = 0
	}
}
