// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package ana

import (
	"math"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/ode"
)

// ColumnFluidPressure computes pressure (p) and intrinsic density (R) of a fluid
// along a column with gravity (g). The numerical solution is:
//
//    R    = R0 + C・(p - p0)   thus   dR/dp = C
//    Z(z) = zmax + T・(z - zmax)   with 0 ≤ T ≤ 1    T is a pseudo variable
//    dZ   = (z - zmax)・dT
//    dp   = R(p)・g・(-dZ)
//    dp   = R(p)・g・(zmax - z)・dT
//    Δz   = zmax - z
//
//            / dp/dT \    / R(p)・g・Δz \
//    dY/dT = |        | = |             |
//            \ dR/dT /    \  C・dp/dT   /
//
type ColumnFluidPressure struct {
	R0  float64 // intrinsic density corresponding to p0
	p0  float64 // pressure corresponding to R0
	C   float64 // compressibility coefficient; e.g. R0/Kbulk or M/(R・θ)
	g   float64 // gravity acceleration (positive constant)
	H   float64 // elevation where (R0,p0) is known
	sol ode.ODE // ODE solver
}

// Init initialises this structure
func (o *ColumnFluidPressure) Init(R0, p0, C, g, H float64, withNum bool) {

	// input data
	o.R0 = R0
	o.p0 = p0
	o.C = C
	o.g = g
	o.H = H

	// numerical solver with ξ := {p, R}
	if withNum {
		silent := true
		o.sol.Init("Radau5", 2, func(f []float64, dT, T float64, ξ []float64, args ...interface{}) error {
			Δz := args[0].(float64)
			R := ξ[1]
			f[0] = R * o.g * Δz // dp/dT
			f[1] = o.C * f[0]   // dR/dT
			return nil
		}, nil, nil, nil, silent)
		o.sol.Distr = false // must be sure to disable this; otherwise it causes problems in parallel runs
	}
}

// Calc computes pressure and density
func (o ColumnFluidPressure) Calc(z float64) (p, R float64) {
	p = o.p0 + (o.R0/o.C)*(math.Exp(o.C*o.g*(o.H-z))-1.0)
	R = o.R0 + o.C*(p-o.p0)
	return
}

// CalcNum computes pressure and density using numerical method
func (o ColumnFluidPressure) CalcNum(z float64) (p, R float64) {
	Δz := o.H - z
	ξ := []float64{o.p0, o.R0}
	err := o.sol.Solve(ξ, 0, 1, 1, false, Δz)
	if err != nil {
		chk.Panic("ColumnFluidPressure failed when calculating pressure using ODE solver: %v", err)
	}
	return ξ[0], ξ[1]
}
