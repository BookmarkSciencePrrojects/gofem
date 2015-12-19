// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mreten

import (
	"math"
	"strings"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
)

// BrooksCorey implements Books and Corey' model
type BrooksCorey struct {

	// parameters
	λ     float64 // slope coefficient
	pcae  float64 // air-entry pressure
	slmin float64 // residual (minimum) saturation
	slmax float64 // maximum saturation
}

// add model to factory
func init() {
	allocators["bc"] = func() Model { return new(BrooksCorey) }
}

// Init initialises model
func (o *BrooksCorey) Init(prms fun.Prms) (err error) {
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
			return chk.Err("bc: parameter named %q is incorrect\n", p.N)
		}
	}
	return
}

// GetPrms gets (an example) of parameters
func (o BrooksCorey) GetPrms(example bool) fun.Prms {
	return []*fun.Prm{
		&fun.Prm{N: "lam", V: 0.5},
		&fun.Prm{N: "pcae", V: 0.2},
		&fun.Prm{N: "slmin", V: 0.1},
		&fun.Prm{N: "slmax", V: 1.0},
	}
}

// SlMin returns sl_min
func (o BrooksCorey) SlMin() float64 {
	return o.slmin
}

// SlMax returns sl_max
func (o BrooksCorey) SlMax() float64 {
	return o.slmax
}

// Sl computes sl directly from pc
func (o BrooksCorey) Sl(pc float64) float64 {
	if pc <= o.pcae {
		return o.slmax
	}
	return o.slmin + (o.slmax-o.slmin)*math.Pow(o.pcae/pc, o.λ)
}

// Cc computes Cc(pc) := dsl/dpc
func (o BrooksCorey) Cc(pc, sl float64, wet bool) (float64, error) {
	if pc <= o.pcae {
		return 0, nil
	}
	return -(o.slmax - o.slmin) * o.λ * math.Pow(o.pcae/pc, o.λ) / pc, nil
}

// L computes L = ∂Cc/∂pc
func (o BrooksCorey) L(pc, sl float64, wet bool) (float64, error) {
	if pc <= o.pcae {
		return 0, nil
	}
	return (o.slmax - o.slmin) * o.λ * (o.λ + 1.0) * math.Pow(o.pcae/pc, o.λ) / (pc * pc), nil
}

// J computes J = ∂Cc/∂sl
func (o BrooksCorey) J(pc, sl float64, wet bool) (float64, error) {
	return 0, nil
}

// Derivs computes all derivatives
func (o BrooksCorey) Derivs(pc, sl float64, wet bool) (L, Lx, J, Jx, Jy float64, err error) {
	if pc <= o.pcae {
		return
	}
	cf := (o.slmax - o.slmin) * o.λ
	pc2 := pc * pc
	pp := math.Pow(o.pcae/pc, o.λ)
	dppdpc := -o.λ * math.Pow(o.pcae/pc, o.λ) / pc
	L = cf * (o.λ + 1.0) * pp / pc2
	Lx = cf * (o.λ + 1.0) * (dppdpc - 2.0*pp/pc) / pc2
	return
}
