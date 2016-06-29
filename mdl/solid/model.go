// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package solid implements models for solids based on continuum mechanics
/*
 *            |    Rate
 *  ============================================
 *            |
 *            | dσdt = f(σ,dεdt)
 *    Small   | σ_(n+1) = σ_(n) + Δt * f_(n+1)
 *            | StressUpdate
 *            | D = dσ/dε_(n+1)
 *            | ConsistentD
 *            |
 *  --------------------------------------------
 *            |
 *    Large   | dσdt = f(σ,F,dFdt)
 *            | D = dσdF_(n+1)
 *            |
 *            |
 */
package solid

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
)

// Model defines the interface for solid models
type Model interface {
	Init(ndim int, pstress bool, prms fun.Prms) error // initialises model
	InitIntVars(σ []float64) (*State, error)          // initialises AND allocates internal (secondary) variables
	GetPrms() fun.Prms                                // gets (an example) of parameters
	GetRho() float64                                  // returns density
	Clean()                                           // clean resources as when calling C code
}

// Small defines rate type solid models for small strain analyses
type Small interface {
	Update(s *State, ε, Δε []float64, eid, ipid int, time float64) error // updates stresses for given strains
	CalcD(D [][]float64, s *State, firstIt bool) error                   // computes D = dσ_new/dε_new consistent with StressUpdate
	ContD(D [][]float64, s *State) error                                 // computes D = dσ_new/dε_new continuous
}

// SmallStrainUpdater define small-strain models that can update strains for given stresses
type SmallStrainUpdater interface {
	StrainUpdate(s *State, Δσ []float64) error // updates strains for given stresses (small strains formulation)
}

// Large defines rate type solid models for large deformation analyses
type Large interface {
	Update(s *State, F, FΔ [][]float64) error              // updates stresses for new deformation F and FΔ
	CalcA(A [][][][]float64, s *State, firstIt bool) error // computes tangent modulus A = (2/J) * ∂τ/∂b . b - σ palm I
}

// OneD specialises Model to 1D
type OneD interface {
	InitIntVars1D() (*OnedState, error)                         // initialises AND allocates internal (secondary) variables
	Update(s *OnedState, ε, Δε, aux float64) error              // update state
	CalcD(s *OnedState, firstIt bool) (float64, float64, error) // computes D = dσ_new/dε_new consistent with StressUpdate
	GetA() float64                                              // returns cross-sectional area
}

// New returns new solid model
func New(name string) (model Model, err error) {
	allocator, ok := allocators[name]
	if !ok {
		return nil, chk.Err("model %q is not available in 'solid' database", name)
	}
	return allocator(), nil
}

// allocators holds all available solid models; modelname => allocator
var allocators = map[string]func() Model{}
