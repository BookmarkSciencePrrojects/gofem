// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package retention

import (
	"math"
	"strings"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
)

// Lin implements a linear retetion model: sl(pc) := 1 - λ*pc
type Lin struct {

	// parameters
	λ     float64 // slope coefficient
	pcae  float64 // air-entry pressure
	slmin float64 // residual (minimum) saturation
	slmax float64 // maximum saturation

	// derived
	pcres float64 // residual pc corresponding to slmin
}

// add model to factory
func init() {
	allocators["lin"] = func() Model { return new(Lin) }
}

// Init initialises model
func (o *Lin) Init(prms fun.Prms) (err error) {
	o.slmax = 1.0
	for _, p := range prms {
		switch strings.ToLower(p.N) {
		case "lam":
			o.λ = p.V
		case "pcae":
			o.pcae = p.V
		case "slmin":
			o.slmin = p.V
		case "slmax":
			o.slmax = p.V
		default:
			return chk.Err("lin: parameter named %q is incorrect\n", p.N)
		}
	}
	if o.λ < 1e-13 {
		o.λ = 0
		o.pcres = math.MaxFloat64
	} else {
		o.pcres = o.pcae + (o.slmax-o.slmin)/o.λ
	}
	return
}

// GetPrms gets (an example) of parameters
func (o Lin) GetPrms(example bool) fun.Prms {
	return []*fun.Prm{
		&fun.Prm{N: "lam", V: 0.5},
		&fun.Prm{N: "pcae", V: 0.2},
		&fun.Prm{N: "slmin", V: 0.1},
		&fun.Prm{N: "slmax", V: 1.0},
	}
}

// SlMin returns sl_min
func (o Lin) SlMin() float64 {
	return o.slmin
}

// SlMax returns sl_max
func (o Lin) SlMax() float64 {
	return o.slmax
}

// Sl computes sl directly from pc
func (o Lin) Sl(pc float64) float64 {
	if pc <= o.pcae {
		return o.slmax
	}
	if pc >= o.pcres {
		return o.slmin
	}
	return o.slmax - o.λ*(pc-o.pcae)
}

// Cc computes Cc(pc) := dsl/dpc
func (o Lin) Cc(pc, sl float64, wet bool) (float64, error) {
	if pc <= o.pcae || pc >= o.pcres {
		return 0, nil
	}
	return -o.λ, nil
}

// L computes L = ∂Cc/∂pc
func (o Lin) L(pc, sl float64, wet bool) (float64, error) {
	return 0, nil
}

// J computes J = ∂Cc/∂sl
func (o Lin) J(pc, sl float64, wet bool) (float64, error) {
	return 0, nil
}

// Derivs compute ∂Cc/∂pc and ∂²Cc/∂pc²
func (o Lin) Derivs(pc, sl float64, wet bool) (L, Lx, J, Jx, Jy float64, err error) {
	return
}
