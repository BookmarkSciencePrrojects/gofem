// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"

	"github.com/cpmech/gofem/inp"
)

// DynCoefs calculates θ-method, Newmark's or HHT coefficients.
//  Notes:
//   θ1  -- Newmark parameter (gamma)  [0 <= θ1 <= 1]
//   θ2  -- Newmark parameter (2*beta) [0 <= θ2 <= 1]
//   HHT -- use Hilber-Hughes-Taylor method ?
//   α   -- Hilber-Hughes-Taylor parameter [-1/3 <= α <= 0]
//   if HHT==True, θ1 and θ2 are automatically calculated for unconditional stability
type DynCoefs struct {

	// input
	θ, θ1, θ2, α float64
	HHT          bool

	// derived
	β1, β2     float64
	α1, α2, α3 float64
	α4, α5, α6 float64
	α7, α8     float64
	hmin       float64
}

// Init initialises this structure
func (o *DynCoefs) Init(dat *inp.SolverData) {

	// hmin
	o.hmin = dat.DtMin

	// HHT
	o.HHT = dat.HHT

	// θ-method
	o.θ = dat.Theta
	if o.θ < 1e-5 || o.θ > 1.0 {
		chk.Panic("θ-method requires 1e-5 <= θ <= 1.0 (θ = %v is incorrect)", o.θ)
	}

	// HHT method
	if dat.HHT {
		o.α = dat.HHTalp
		if o.α < -1.0/3.0 || o.α > 0.0 {
			chk.Panic("HHT method requires: -1/3 <= α <= 0 (α = %v is incorrect)", o.α)
		}
		o.θ1 = (1.0 - 2.0*o.α) / 2.0
		o.θ2 = (1.0 - o.α) * (1.0 - o.α) / 2.0

		// Newmark's method
	} else {
		o.θ1, o.θ2 = dat.Theta1, dat.Theta2
		if o.θ1 < 0.0001 || o.θ1 > 1.0 {
			chk.Panic("θ1 must be between 0.0001 and 1.0 (θ1 = %v is incorrect)", o.θ1)
		}
		if o.θ2 < 0.0001 || o.θ2 > 1.0 {
			chk.Panic("θ2 must be between 0.0001 and 1.0 (θ2 = %v is incorrect)", o.θ2)
		}
	}
}

// CalcBoth computes betas and alphas
func (o *DynCoefs) CalcBoth(Δt float64) (err error) {
	err = o.CalcBetas(Δt)
	if err != nil {
		return
	}
	err = o.CalcAlphas(Δt)
	return
}

// CalcBetas computes only betas
func (o *DynCoefs) CalcBetas(Δt float64) (err error) {

	// timestep
	h := Δt
	if h < o.hmin {
		return chk.Err("θ-method requires h >= %v (h = %v is incorrect)", o.hmin, h)
	}

	// β coefficients
	o.β1 = 1.0 / (o.θ * h)
	o.β2 = (1.0 - o.θ) / o.θ
	return
}

// CalcAlphas computes only alphas
func (o *DynCoefs) CalcAlphas(Δt float64) (err error) {

	// timestep
	h := Δt
	if h < o.hmin {
		return chk.Err("Newmark/HHT method requires h >= %v (h = %v is incorrect)", o.hmin, h)
	}

	// α coefficients
	H := h * h / 2.0
	o.α1, o.α2, o.α3 = 1.0/(o.θ2*H), h/(o.θ2*H), 1.0/o.θ2-1.0
	o.α4, o.α5, o.α6 = o.θ1*h/(o.θ2*H), 2.0*o.θ1/o.θ2-1.0, (o.θ1/o.θ2-1.0)*h

	// HHT method
	o.α7 = o.α4
	o.α8 = 1.0
	if o.HHT {
		o.α7, o.α8 = (1.0+o.α)*o.α4, 1.0+o.α
	}
	return
}

// Print prints coefficients
func (o *DynCoefs) Print() {
	io.Pfgrey("θ=%v, θ1=%v, θ2=%v, α=%v\n", o.θ, o.θ1, o.θ2, o.α)
	io.Pfgrey("HHT=%v\n", o.HHT)
	io.Pfgrey("β1=%v, β2=%v\n", o.β1, o.β2)
	io.Pfgrey("α1=%v, α2=%v, α3=%v, α4=%v, α5=%v, α6=%v\n", o.α1, o.α2, o.α3, o.α4, o.α5, o.α6)
	io.Pfgrey("α7=%v, α8=%v\n", o.α7, o.α8)
}
