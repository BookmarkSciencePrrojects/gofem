// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package retention implements models for liquid retention curves
//  References:
//   [1] Pedroso DM, Sheng D and Zhao, J (2009) The concept of reference curves for constitutive
//       modelling in soil mechanics, Computers and Geotechnics, 36(1-2), 149-165,
//       http://dx.doi.org/10.1016/j.compgeo.2008.01.009
//   [2] Pedroso DM and Williams DJ (2010) A novel approach for modelling soil-water
//       characteristic curves with hysteresis, Computers and Geotechnics, 37(3), 374-380,
//       http://dx.doi.org/10.1016/j.compgeo.2009.12.004
//   [3] Pedroso DM and Williams DJ (2011) Automatic Calibration of soil-water characteristic
//       curves using genetic algorithms. Computers and Geotechnics, 38(3), 330-340,
//       http://dx.doi.org/10.1016/j.compgeo.2010.12.004
package retention

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/ode"
)

// Model implements a liquid retention model (LRM)
//  Derivs computes (see [1] page 618):
//    L  = ∂Cc/∂pc
//    Lx = ∂²Cc/∂pc²
//    J  = ∂Cc/∂sl
//    Jx == ∂²Cc/(∂pc ∂sl)
//    Jy == ∂²Cc/∂sl²
//  References:
//   [1] Pedroso DM (2015) A consistent u-p formulation for porous media with hysteresis.
//       Int Journal for Numerical Methods in Engineering, 101(8) 606-634
//       http://dx.doi.org/10.1002/nme.4808
type Model interface {
	Init(prms fun.Prms) error                                              // initialises retention model
	GetPrms(example bool) fun.Prms                                         // gets (an example) of parameters
	SlMin() float64                                                        // returns sl_min
	SlMax() float64                                                        // returns sl_max
	Cc(pc, sl float64, wet bool) (float64, error)                          // computes Cc = f = ∂sl/∂pc
	L(pc, sl float64, wet bool) (float64, error)                           // computes L = ∂Cc/∂pc
	J(pc, sl float64, wet bool) (float64, error)                           // computes J = ∂Cc/∂sl
	Derivs(pc, sl float64, wet bool) (L, Lx, J, Jx, Jy float64, err error) // computes all derivatives
}

// Nonrate is a subset of LRM that directly computes saturation from capillary pressure
type Nonrate interface {
	Sl(pc float64) float64 // compute sl directly from pc
}

// Update updates pc and sl for given Δpc. An implicit ODE solver is used.
func Update(mdl Model, pc0, sl0, Δpc float64) (slNew float64, err error) {

	// wetting flag
	wet := Δpc < 0

	// callback functions
	//   x      = [0.0, 1.0]
	//   pc     = pc0 + x * Δpc
	//   y[0]   = sl
	//   f(x,y) = dy/dx = dsl/dpc * dpc/dx = Cc * Δpc
	//   J(x,y) = df/dy = DCcDsl * Δpc
	fcn := func(f []float64, dx, x float64, y []float64) (e error) {
		f[0], e = mdl.Cc(pc0+x*Δpc, y[0], wet)
		f[0] *= Δpc
		return nil
	}
	jac := func(dfdy *la.Triplet, dx, x float64, y []float64) (e error) {
		if dfdy.Max() == 0 {
			dfdy.Init(1, 1, 1)
		}
		J, e := mdl.J(pc0+x*Δpc, y[0], wet)
		dfdy.Start()
		dfdy.Put(0, 0, J)
		return
	}

	// ode solver
	var odesol ode.Solver
	odesol.Init("Radau5", 1, fcn, jac, nil, nil)
	odesol.SetTol(1e-10, 1e-7)
	odesol.Distr = false // this is important to avoid problems with MPI runs

	// solve
	y := []float64{sl0}
	err = odesol.Solve(y, 0, 1, 1, false)
	slNew = y[0]
	return
}

// New returns new liquid retention model
func New(name string) (model Model, err error) {
	allocator, ok := allocators[name]
	if !ok {
		return nil, chk.Err("model %q is not available in 'retention' database", name)
	}
	return allocator(), nil
}

// allocators holds all available models
var allocators = map[string]func() Model{}
