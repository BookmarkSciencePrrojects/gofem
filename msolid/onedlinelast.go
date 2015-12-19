// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import "github.com/cpmech/gosl/fun"

// OnedLinElast implements a linear elastic model for 1D elements
type OnedLinElast struct {
	E   float64 // Young's modulus
	G   float64 // shear modulus
	A   float64 // cross-sectional area
	I22 float64 // moment of inertia of cross section about y2-axis
	I11 float64 // moment of inertia of cross section about y1-axis
	Jtt float64 // torsional constant
	Rho float64 // density
}

// add model to factory
func init() {
	allocators["oned-elast"] = func() Model { return new(OnedLinElast) }
}

// Clean clean resources
func (o *OnedLinElast) Clean() {
}

// Init initialises model
func (o *OnedLinElast) Init(ndim int, pstress bool, prms fun.Prms) (err error) {
	for _, p := range prms {
		switch p.N {
		case "E":
			o.E = p.V
		case "G":
			o.G = p.V
		case "A":
			o.A = p.V
		case "I22":
			o.I22 = p.V
		case "I11":
			o.I11 = p.V
		case "Jtt":
			o.Jtt = p.V
		case "rho":
			o.Rho = p.V
		}
	}
	return
}

// GetPrms gets (an example) of parameters
func (o OnedLinElast) GetPrms() fun.Prms {
	return []*fun.Prm{
		&fun.Prm{N: "E", V: 2.0000e+08},
		&fun.Prm{N: "G", V: 7.5758e+07},
		&fun.Prm{N: "A", V: 1.0000e-02},
		&fun.Prm{N: "I22", V: 8.3333e-06},
		&fun.Prm{N: "I11", V: 8.3333e-06},
		&fun.Prm{N: "Jtt", V: 1.4063e-05},
		&fun.Prm{N: "rho", V: 7.8500e+00},
	}
}

// InitIntVars initialises internal (secondary) variables
func (o OnedLinElast) InitIntVars1D() (s *OnedState, err error) {
	s = NewOnedState(0, 0)
	return
}

// Update updates stresses for given strains
func (o OnedLinElast) Update(s *OnedState, ε, Δε float64) (err error) {
	s.Sig += o.E * Δε
	return
}

// CalcD computes D = dσ_new/dε_new consistent with StressUpdate
func (o OnedLinElast) CalcD(s *OnedState, firstIt bool) (float64, error) {
	return o.E, nil
}
