// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package diffusion

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// M1 implements a model for diffusion problems with nonlinear coefficient
//
//   kten = kval(u) * kcte
//
//   kval = a0  +  a1 u  +  a2 u² +  a3 u³
//
type M1 struct {
	a0, a1, a2, a3 float64
	Rho            float64
	Kcte           [][]float64
}

// add model to factory
func init() {
	allocators["m1"] = func() Model { return new(M1) }
}

// Init initialises this structure
func (o *M1) Init(ndim int, prms fun.Prms) (err error) {

	// a[i] parameters
	prms.Connect(&o.a0, "a0", "a0 M1 model")
	prms.Connect(&o.a1, "a1", "a1 M1 model")
	prms.Connect(&o.a2, "a2", "a2 M1 model")
	prms.Connect(&o.a3, "a3", "a3 M1 model")
	prms.Connect(&o.Rho, "rho", "rho M1 model")

	// keys
	keys := []string{"kx", "ky"}
	if ndim == 3 {
		keys = []string{"kx", "ky", "kz"}
	}

	// kcte parameters
	var kx, ky, kz float64
	k_values, k_found := prms.GetValues(keys)
	if !utl.BoolAllTrue(k_found) {
		p := prms.Find("k")
		if p == nil {
			return chk.Err("M1 model: either 'k' (isotropic) or ['kx', 'ky', 'kz'] must be given in database of material parameters")
		}
		kx, ky, kz = p.V, p.V, p.V
	} else {
		kx, ky = k_values[0], k_values[1]
		if ndim == 3 {
			kz = k_values[2]
		}
	}

	// ktensor
	o.Kcte = la.MatAlloc(ndim, ndim)
	o.Kcte[0][0] = kx
	o.Kcte[1][1] = ky
	if ndim == 3 {
		o.Kcte[2][2] = kz
	}
	return
}

// Kval computes k(u)
func (o *M1) Kval(u float64) float64 {
	return o.a0 + o.a1*u + o.a2*u*u + o.a3*u*u*u
}

// DkDu computes dk/du
func (o *M1) DkDu(u float64) float64 {
	return o.a1 + 2.0*o.a2*u + 3.0*o.a3*u*u
}

// Kten computes ktensor = kval(u) * kcte
func (o *M1) Kten(kten [][]float64, u float64) {
	la.MatCopy(kten, o.Kval(u), o.Kcte)
}
