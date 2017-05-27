// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package fluid implements models for fluid density
package fluid

import (
	"math"

	"github.com/cpmech/gosl/fun"
)

// Model implements a model to compute pressure (p) and intrinsic density (R) of a fluid
// along a column with gravity (g). The model is:
//   R(p) = R0 + C・(p - p0)   thus   dR/dp = C
type Model struct {

	// material data
	R0  float64 // intrinsic density corresponding to p0
	P0  float64 // pressure corresponding to R0
	C   float64 // compressibility coefficient; e.g. R0/Kbulk or M/(R・θ)
	Gas bool    // is gas instead of liquid?

	// additional data
	H    float64 // elevation where (R0,p0) is known
	Grav float64 // gravity acceleration (positive constant)
}

// Init initialises this structure
func (o *Model) Init(prms fun.Params, H, grav float64) {
	for _, p := range prms {
		switch p.N {
		case "R0":
			o.R0 = p.V
		case "P0":
			o.P0 = p.V
		case "C":
			o.C = p.V
		case "gas":
			o.Gas = p.V > 0
		}
	}
	o.H = H
	o.Grav = grav
}

// GetPrms gets (an example of) parameters
//  Input:
//   example -- returns example of parameters; othewise returs current parameters
//  Note:
//   Gas variable is used to return dry air properties instead of water
func (o Model) GetPrms(example bool) fun.Params {
	if example {
		if o.Gas {
			return fun.Params{ // dry air
				&fun.P{N: "R0", V: 0.0012}, // [Mg/m³]
				&fun.P{N: "P0", V: 0.0},    // [kPa]
				&fun.P{N: "C", V: 1.17e-5}, // [Mg/(m³・kPa)]
				&fun.P{N: "Gas", V: 1},     // [-]
			}
		}
		return fun.Params{ // water
			&fun.P{N: "R0", V: 1.0},    // [Mg/m³]
			&fun.P{N: "P0", V: 0.0},    // [kPa]
			&fun.P{N: "C", V: 4.53e-7}, // [Mg/(m³・kPa)]
			&fun.P{N: "Gas", V: 0},     // [-]
		}
	}
	var gas float64
	if o.Gas {
		gas = 1
	}
	return fun.Params{
		&fun.P{N: "R0", V: o.R0},
		&fun.P{N: "P0", V: o.P0},
		&fun.P{N: "C", V: o.C},
		&fun.P{N: "Gas", V: gas},
	}
}

// Calc computes pressure and density
func (o Model) Calc(z float64) (p, R float64) {
	p = o.P0 + (o.R0/o.C)*(math.Exp(o.C*o.Grav*(o.H-z))-1.0)
	R = o.R0 + o.C*(p-o.P0)
	return
}

// Plot plots pressure and density along height of column
func (o Model) Plot(dirout, fnkey string, np int) {

	/* TODO
	Z := utl.LinSpace(0, o.H, np)
	P := make([]float64, np)
	R := make([]float64, np)
	for i, z := range Z {
		P[i], R[i] = o.Calc(z)
	}

	pMaxLin := o.R0 * o.Grav * o.H
	subscript := "\\ell"
	if o.Gas {
		subscript = "g"
	}

	plt.Subplot(2, 1, 1)
	plt.Plot(P, Z, "'b-', clip_on=0")
	plt.Plot([]float64{o.P0, pMaxLin}, []float64{o.H, 0}, "'k--', color='gray'")
	plt.Gll("$p_{"+subscript+"}$", "$z$", "")

	plt.Subplot(2, 1, 2)
	plt.Plot(R, Z, "'r-', clip_on=0")
	plt.Plot([]float64{o.R0, o.R0 + o.C*pMaxLin}, []float64{o.H, 0}, "'k--', color='gray'")
	plt.Gll("$\\rho_{"+subscript+"}$", "$z$", "")
	plt.SetTicksNormal()

	plt.SaveD(dirout, fnkey+".eps")
	*/
}
