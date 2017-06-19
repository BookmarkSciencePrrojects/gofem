// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package ele implements finite elements
package ele

import (
	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/fun/dbf"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// Element defines what all elements must implement
type Element interface {

	// information and initialisation
	Id() int                                             // returns the cell Id
	SetEqs(eqs [][]int, mixedform_eqs []int) (err error) // set equations

	// conditions (natural BCs and element's)
	SetEleConds(key string, f dbf.T, extra string) (err error) // set element conditions

	// called for each time step
	InterpStarVars(sol *Solution) (err error) // interpolate star variables to integration points

	// called for each iteration
	AddToRhs(fb []float64, sol *Solution) (err error)                // adds -R to global residual vector fb
	AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (err error) // adds element K to global Jacobian matrix Kb

	// reading and writing of element data
	Encode(enc utl.Encoder) (err error) // encodes internal variables
	Decode(dec utl.Decoder) (err error) // decodes internal variables
}

// WithIntVars defines elements with {z,q} internal variables
type WithIntVars interface {
	Update(sol *Solution) (err error)                              // perform (tangent) update
	SetIniIvs(sol *Solution, ivs map[string][]float64) (err error) // sets initial ivs for given values in sol and ivs map
	BackupIvs(aux bool) (err error)                                // create copy of internal variables
	RestoreIvs(aux bool) (err error)                               // restore internal variables from copies
	Ureset(sol *Solution) (err error)                              // fixes internal variables after u (displacements) have been zeroed
}

// Connector defines connector elements; elements that depend upon others
type Connector interface {
	Id() int                                                       // returns the cell Id
	Connect(cid2elem []Element, c *inp.Cell) (nnzK int, err error) // connect multiple elements; e.g.: connect rod/solid elements in Rjoints
}

// CanExtrapolate defines elements with functions to extrapolate internal values
type CanExtrapolate interface {
	AddToExt(sol *Solution) (err error) // adds extrapolated values to global array
}

// CanOutputIps defines elements that can output integration points' values
type CanOutputIps interface {
	Id() int                            // returns the cell Id
	OutIpCoords() [][]float64           // coordinates of integration points
	OutIpKeys() []string                // integration points' keys; e.g. "pl", "sl"
	OutIpVals(M *IpsMap, sol *Solution) // integration points' values corresponding to keys
}

// WithFixedKM defines elements with fixed K,M matrices; to be recomputed if prms are changed
type WithFixedKM interface {
	Recompute(withM bool) // recompute K and M
}
