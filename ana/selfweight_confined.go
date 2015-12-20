// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package ana implements analytical solutions
package ana

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
)

// ConfinedSelfWeight computes the solution to a simple confined linear elastic domain under gravity
//
//     ▷ o-----------o ◁
//     ▷ |           | ◁
//     ▷ |    E, ρ   | ◁       negative stress means compression
//  h  ▷ |    ν, g   | ◁       g = 10  =>  b = -10 * ρ (body force)
//     ▷ |           | ◁
//     ▷ o-----------o ◁
//       △  △  △  △  △
//             w
type ConfinedSelfWeight struct {
	// input
	E float64 // Young's modulus
	ν float64 // Poisson's coefficient
	ρ float64 // density
	g float64 // gravity constant (positive value)
	h float64 // height
	w float64 // width

	// derived
	d float64 // auxiliary coefficient = ν/(1-ν)   = E/((1+ν)(1-2ν))
	M float64 // P-wave modulus
}

// Init initialises this structure
func (o *ConfinedSelfWeight) Init(prms fun.Prms) {

	// default values
	o.E = 1000.0
	o.ν = 0.25
	o.ρ = 2.0
	o.g = 10.0
	o.h = 1.0
	o.w = 1.0

	// parameters
	for _, p := range prms {
		switch p.N {
		case "E":
			o.E = p.V
		case "nu":
			o.ν = p.V
		case "rho":
			o.ρ = p.V
		case "g":
			o.g = p.V
		case "h":
			o.h = p.V
		case "w":
			o.w = p.V
		}
	}

	// derived
	//o.c = o.E / ((1.0 + o.ν) * (1.0 - 2.0*o.ν))
	//o.M = o.c * (1.0 - o.ν)
	o.d = o.ν / (1.0 - o.ν)
	o.M = o.E * (1.0 - o.ν) / ((1.0 + o.ν) * (1.0 - 2.0*o.ν))
}

// Stress computes stress components
func (o ConfinedSelfWeight) Stress(t float64, x []float64) (σ []float64) {
	ndim := len(x)             // space dimension
	z := x[ndim-1]             // elevation
	b := o.g * t               // body force
	σv := -o.ρ * b * (o.h - z) // vertical stress
	σh := o.d * σv             // horizontal stress
	σ = make([]float64, 2*ndim)
	if ndim == 2 {
		σ[0], σ[1], σ[2] = σh, σv, σh
		return
	}
	σ[0], σ[1], σ[2] = σh, σh, σv
	return
}

// Displ computes displacement components
func (o ConfinedSelfWeight) Displ(t float64, x []float64) (u []float64) {
	ndim := len(x)              // space dimension
	z := x[ndim-1]              // elevation
	b := o.g * t                // body force
	α := -o.ρ * b / o.M         // auxiliary coefficient
	uv := α * (o.h - z/2.0) * z // vertical displacement
	u = make([]float64, ndim)
	if ndim == 2 {
		u[0], u[1], u[2] = 0, uv, 0
		return
	}
	u[0], u[1], u[2] = 0, 0, uv
	return
}

// CheckStress check stresses
func (o ConfinedSelfWeight) CheckStress(tst *testing.T, t float64, σ, x []float64, tol float64) {
	σana := o.Stress(t, x)
	chk.Vector(tst, "σ", tol, σ, σana)
}

// CheckDispl checks displacements
func (o ConfinedSelfWeight) CheckDispl(tst *testing.T, t float64, u, x []float64, tol float64) {
	uana := o.Displ(t, x)
	chk.Vector(tst, "u", tol, u, uana)
}
