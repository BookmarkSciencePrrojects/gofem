// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package msolid implements models for solids based on continuum mechanics
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
package msolid

import (
	"log"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
)

// Model defines the interface for solid models
type Model interface {
	Init(ndim int, pstress bool, prms fun.Prms) error // initialises model
	GetPrms() fun.Prms                                // gets (an example) of parameters
	InitIntVars(σ []float64) (*State, error)          // initialises AND allocates internal (secondary) variables
}

// Small defines rate type solid models for small strain analyses
type Small interface {
	Update(s *State, ε, Δε []float64, eid, ipid int, time float64) error // updates stresses for given strains
	CalcD(D [][]float64, s *State, firstIt bool) error                   // computes D = dσ_new/dε_new consistent with StressUpdate
	ContD(D [][]float64, s *State) error                                 // computes D = dσ_new/dε_new continuous
}

// Large defines rate type solid models for large deformation analyses
type Large interface {
	Update(s *State, F, FΔ [][]float64) error              // updates stresses for new deformation F and FΔ
	CalcA(A [][][][]float64, s *State, firstIt bool) error // computes tangent modulus A = (2/J) * ∂τ/∂b . b - σ palm I
}

// SmallStrainUpdater define small-strain models that can update strains for given stresses
type SmallStrainUpdater interface {
	StrainUpdate(s *State, Δσ []float64) error // updates strains for given stresses (small strains formulation)
}

// GetModel returns (existent or new) solid model
//  simfnk    -- unique simulation filename key
//  matname   -- name of material
//  modelname -- model name
//  getnew    -- force a new allocation; i.e. do not use any model found in database
//  Note: returns model==nil on errors
func GetModel(simfnk, matname, modelname string, getnew bool) (model Model, existent bool) {

	// get new model, regardless wheter it exists in database or not
	if getnew {
		allocator, ok := allocators[modelname]
		if !ok {
			return nil, false
		}
		return allocator(), false
	}

	// search database
	key := io.Sf("%s_%s_%s", simfnk, matname, modelname)
	if model, ok := _models[key]; ok {
		return model, true
	}

	// if not found, get new
	allocator, ok := allocators[modelname]
	if !ok {
		return nil, false
	}
	model = allocator()
	_models[key] = model
	return model, false
}

// LogModels prints to log information on existent and allocated Models
func LogModels() {
	l := "msolid: available:"
	for name, _ := range allocators {
		l += " " + name
	}
	log.Println(l)
	l = "msolid: allocated:"
	for key, _ := range _models {
		l += " " + io.Sf("%q", key)
	}
	log.Println(l)
	onedLogModels()
}

// allocators holds all available solid models; modelname => allocator
var allocators = map[string]func() Model{}

// _models holds pre-allocated solid models (internal); key => Solid
var _models = map[string]Model{}
