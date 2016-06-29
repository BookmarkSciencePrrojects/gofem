// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package solid

import "github.com/cpmech/gosl/fun"

// LinElast implements a linear elastic model
type LinElast struct {
	SmallElasticity
}

// add model to factory
func init() {
	allocators["lin-elast"] = func() Model { return new(LinElast) }
}

// Clean clean resources
func (o *LinElast) Clean() {
}

// Init initialises model
func (o *LinElast) Init(ndim int, pstress bool, prms fun.Prms) (err error) {
	return o.SmallElasticity.Init(ndim, pstress, prms)
}

// GetPrms gets (an example) of parameters
func (o LinElast) GetPrms() fun.Prms {
	return o.SmallElasticity.GetPrms()
}

// InitIntVars initialises internal (secondary) variables
func (o LinElast) InitIntVars(σ []float64) (s *State, err error) {
	s = NewState(o.Nsig, 0, false, false)
	copy(s.Sig, σ)
	return
}

// Update updates stresses for given strains
func (o LinElast) Update(s *State, ε, Δε []float64, eid, ipid int, time float64) (err error) {
	return o.SmallElasticity.Update(s, Δε)
}

// CalcD computes D = dσ_new/dε_new consistent with StressUpdate
func (o LinElast) CalcD(D [][]float64, s *State, firstIt bool) (err error) {
	return o.SmallElasticity.CalcD(D, s)
}

// ContD computes D = dσ_new/dε_new continuous
func (o LinElast) ContD(D [][]float64, s *State) (err error) {
	return o.SmallElasticity.CalcD(D, s)
}
