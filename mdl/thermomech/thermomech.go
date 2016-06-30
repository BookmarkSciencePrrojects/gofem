// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package thermomech

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// Thermomech implements a model for thermomechanical problems with nonlinear conduction coefficient
//
//   kten = kval(u) * kcte
//
//   kval = a0  +  a1 u  +  a2 u² +  a3 u³
//
type Thermomech struct {
	a0, a1, a2, a3, H 	float64
	Bcte			[]float64
	Kcte			[][]float64
}

// add model to factory
func init() {
	allocators["thermomech"] = func() Model { return new(Thermomech) }
}

// Init initialises this structure
func (o *Thermomech) Init(ndim int, prms fun.Prms) (err error) {

	// a[i] parameters
	prms.Connect(&o.a0, "a0", "a0 Thermomech model")
	prms.Connect(&o.a1, "a1", "a1 Thermomech model")
	prms.Connect(&o.a2, "a2", "a2 Thermomech model")
	prms.Connect(&o.a3, "a3", "a3 Thermomech model")
	prms.Connect(&o.H, "h", "h Thermomech model")

	// expansion keys
	b_keys := []string{"bx", "by"}
	if ndim == 3 {
		b_keys = []string{"bx", "by", "bz"}
	}

	// bcte parameters
	var bx, by, bz float64
	b_values, b_found := prms.GetValues(b_keys)
	if !utl.BoolAllTrue(b_found) {
		p := prms.Find("b")
		if p == nil {
			return chk.Err("Thermomech model: either 'b' (isotropic) or ['bx', 'by', 'bz'] must be given in database of material parameters")
		}
		bx, by, bz = p.V, p.V, p.V
	} else {
		bx, by = b_values[0], b_values[1]
		if ndim == 3 {
			bz = b_values[2]
		}
	}

	// btensor
	o.Bcte = make([]float64, 3)
	o.Bcte[0] = bx
	o.Bcte[1] = by
	if ndim == 3 {
		o.Bcte[2] = bz
	}

	// conductivity keys
	k_keys := []string{"kx", "ky"}
	if ndim == 3 {
		k_keys = []string{"kx", "ky", "kz"}
	}

	// kcte parameters
	var kx, ky, kz float64
	k_values, k_found := prms.GetValues(k_keys)
	if !utl.BoolAllTrue(k_found) {
		p := prms.Find("k")
		if p == nil {
			return chk.Err("Thermomech model: either 'k' (isotropic) or ['kx', 'ky', 'kz'] must be given in database of material parameters")
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
func (o *Thermomech) Kval(u float64) float64 {
	return o.a0 + o.a1*u + o.a2*u*u + o.a3*u*u*u
}

// DkDu computes dk/du
func (o *Thermomech) DkDu(u float64) float64 {
	return o.a1 + 2.0*o.a2*u + 3.0*o.a3*u*u
}

// Kten computes ktensor = kval(u) * kcte
func (o *Thermomech) Kten(kten [][]float64, u float64) {
	la.MatCopy(kten, o.Kval(u), o.Kcte)
}
