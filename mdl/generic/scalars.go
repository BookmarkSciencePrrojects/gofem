// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package generic

import "github.com/cpmech/gosl/fun"

// Scalars is a placeholder model to collect scalar parameters
type Scalars struct {
	A, B, C, D, E, F, G, H    float64
	I, J, K, L, M, N, O, P, Q float64
	R, S, T, U, V, W, X, Y, Z float64
}

// add model to factory
func init() {
	allocators["scalars"] = func() Model { return new(Scalars) }
}

// Init initialises this structure
func (o *Scalars) Init(ndim int, prms fun.Params) (err error) {

	// parameters
	prms.Connect(&o.A, "A", "A parameter. Scalars model")
	prms.Connect(&o.B, "B", "B parameter. Scalars model")
	prms.Connect(&o.C, "C", "C parameter. Scalars model")
	prms.Connect(&o.D, "D", "D parameter. Scalars model")
	prms.Connect(&o.E, "E", "E parameter. Scalars model")
	prms.Connect(&o.F, "F", "F parameter. Scalars model")
	prms.Connect(&o.G, "G", "G parameter. Scalars model")
	prms.Connect(&o.H, "H", "H parameter. Scalars model")
	prms.Connect(&o.I, "I", "I parameter. Scalars model")
	prms.Connect(&o.J, "J", "J parameter. Scalars model")
	prms.Connect(&o.K, "K", "K parameter. Scalars model")
	prms.Connect(&o.L, "L", "L parameter. Scalars model")
	prms.Connect(&o.M, "M", "M parameter. Scalars model")
	prms.Connect(&o.N, "N", "N parameter. Scalars model")
	prms.Connect(&o.O, "O", "O parameter. Scalars model")
	prms.Connect(&o.P, "P", "P parameter. Scalars model")
	prms.Connect(&o.Q, "Q", "Q parameter. Scalars model")
	prms.Connect(&o.R, "R", "R parameter. Scalars model")
	prms.Connect(&o.S, "S", "S parameter. Scalars model")
	prms.Connect(&o.T, "T", "T parameter. Scalars model")
	prms.Connect(&o.U, "U", "U parameter. Scalars model")
	prms.Connect(&o.V, "V", "V parameter. Scalars model")
	prms.Connect(&o.W, "W", "W parameter. Scalars model")
	prms.Connect(&o.X, "X", "X parameter. Scalars model")
	prms.Connect(&o.Y, "Y", "Y parameter. Scalars model")
	prms.Connect(&o.Z, "Z", "Z parameter. Scalars model")
	return
}
