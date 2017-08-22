// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package solid

import (
	"math"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/fun/dbf"
)

// RjointM1 implements a 1D plasticity model for rod-joints (links/interface)
//  Note: σc has opposite sign convention: positive means compressive
type RjointM1 struct {
	A_ks  float64 // elasticity constant
	A_τy0 float64 // initial yield stress
	A_kh  float64 // hardening modulus
	A_μ   float64 // friction coefficient
	A_h   float64 // perimeter of beam element
	A_kl  float64 // lateral stiffness
}

// add model to factory
func init() {
	allocators["rjoint-m1"] = func() Model { return new(RjointM1) }
}

// Set_mu sets μ parameter
func (o *RjointM1) Set_mu(mu float64) {
	o.A_μ = mu
}

// Free frees memory
func (o *RjointM1) Free() {
}

// GetRho returns density
func (o *RjointM1) GetRho() float64 {
	return 0
}

// Init initialises model
func (o *RjointM1) Init(ndim int, pstress bool, prms dbf.Params) (err error) {
	for _, p := range prms {
		switch p.N {
		case "ks":
			o.A_ks = p.V
		case "tauy0":
			o.A_τy0 = p.V
		case "kh":
			o.A_kh = p.V
		case "mu":
			o.A_μ = p.V
		case "h":
			o.A_h = p.V
		case "kl":
			o.A_kl = p.V
		}
	}
	ZERO := 1e-7
	if o.A_ks < ZERO || o.A_τy0 < ZERO || o.A_μ < ZERO || o.A_h < ZERO || o.A_kl < ZERO {
		return chk.Err("invalid parameters: {ks=%g, tauy0=%g, mu=%g, h=%g, kl=%g} must be all > 0", o.A_ks, o.A_τy0, o.A_μ, o.A_h, o.A_kl)
	}
	return
}

// GetPrms gets (an example) of parameters
func (o RjointM1) GetPrms() dbf.Params {
	return []*fun.P{
		&fun.P{N: "ks", V: 1e4},
		&fun.P{N: "tauy0", V: 20},
		&fun.P{N: "kh", V: 0},
		&fun.P{N: "mu", V: 0.5},
		&fun.P{N: "h", V: 0.1},
		&fun.P{N: "kl", V: 1e4},
	}
}

// InitIntVars: unused
func (o *RjointM1) InitIntVars(σ []float64) (s *State, err error) {
	return
}

// InitIntVars initialises internal (secondary) variables
func (o RjointM1) InitIntVars1D() (s *OnedState, err error) {
	s = NewOnedState(1, 2) // 1:{ωpb}  2:{q1,q2}
	return
}

// Update updates stresses for given strains
//  Note: σc has opposite sign convention: positive means compressive
func (o *RjointM1) Update(s *OnedState, σcNew, Δω float64) (err error) {

	// limit σcNew
	if σcNew < 0 {
		σcNew = 0
	}

	// internal values
	τ := &s.Sig
	ωpb := &s.Alp[0]

	// trial stress
	τ_tr := (*τ) + o.A_ks*Δω
	f_tr := math.Abs(τ_tr) - (o.A_τy0 + o.A_kh*(*ωpb) + o.A_μ*σcNew)

	// elastic update
	if f_tr <= 0.0 {
		*τ = τ_tr
		s.Loading = false
		return
	}

	// plastic update
	Δγ := f_tr / (o.A_ks + o.A_kh)
	*τ = τ_tr - o.A_ks*Δγ*fun.Sign(τ_tr)
	*ωpb += Δγ
	s.Loading = true
	return
}

// CalcD computes D = dσ_new/dε_new consistent with StressUpdate
func (o *RjointM1) CalcD(s *OnedState, firstIt bool) (DτDω, DτDσc float64, err error) {

	// elastic
	if !s.Loading {
		return o.A_ks, 0, nil
	}

	// plastic
	τ := s.Sig
	DτDω = o.A_ks * o.A_kh / (o.A_ks + o.A_kh)
	DτDσc = o.A_ks * o.A_μ * fun.Sign(τ) / (o.A_ks + o.A_kh)
	return
}
