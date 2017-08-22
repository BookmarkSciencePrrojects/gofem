// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package retention

import (
	"math"
	"strings"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/fun/dbf"
)

// VanGen implements van Genuchten's model
type VanGen struct {

	// parameters
	α, m, n float64 // parameters
	slmin   float64 // minimum sl
	slmax   float64 // maximum sl
	pcmin   float64 // pc limit to consider zero slope
	pclim   float64 // pc limit corresponding to slmin
}

// add model to factory
func init() {
	allocators["vg"] = func() Model { return new(VanGen) }
}

// Init initialises model
func (o *VanGen) Init(prms dbf.Params) (err error) {
	o.pcmin, o.slmax = 1e-3, 1.0
	for _, p := range prms {
		switch strings.ToLower(p.N) {
		case "alp":
			o.α = p.V
		case "m":
			o.m = p.V
		case "n":
			o.n = p.V
		case "slmin":
			o.slmin = p.V
		case "slmax":
			o.slmax = p.V
		case "pcmin":
			o.pcmin = p.V
		default:
			return chk.Err("vg: parameter named %q is incorrect\n", p.N)
		}
	}
	if o.slmin > 0 {
		k := (o.slmax - o.slmin) / o.slmin
		o.pclim = math.Pow((math.Pow(k, 1.0/o.m)-1.0)/math.Pow(o.α, o.n), 1.0/o.n)
	} else {
		o.pclim = 1e+30
	}
	return
}

// GetPrms gets (an example) of parameters
func (o VanGen) GetPrms(example bool) dbf.Params {
	return []*fun.P{
		&fun.P{N: "alp", V: 0.08},
		&fun.P{N: "m", V: 4},
		&fun.P{N: "n", V: 4},
		&fun.P{N: "slmin", V: 0.01},
		&fun.P{N: "slmax", V: 1.0},
		&fun.P{N: "pcmin", V: 1e-3},
	}
}

// SlMin returns sl_min
func (o VanGen) SlMin() float64 {
	return o.slmin
}

// SlMax returns sl_max
func (o VanGen) SlMax() float64 {
	return o.slmax
}

// Sl computes sl directly from pc
func (o VanGen) Sl(pc float64) float64 {
	if pc <= o.pcmin {
		return o.slmax
	}
	if pc >= o.pclim {
		return o.slmin
	}
	c := math.Pow(o.α*pc, o.n)
	fac := o.slmax - o.slmin
	return fac * math.Pow(1+c, -o.m)
}

// Cc computes Cc(pc) := dsl/dpc
func (o VanGen) Cc(pc, sl float64, wet bool) (float64, error) {
	if pc <= o.pcmin || pc >= o.pclim {
		return 0, nil
	}
	c := math.Pow(o.α*pc, o.n)
	fac := o.slmax - o.slmin
	return -fac * c * math.Pow(c+1.0, -o.m-1.0) * o.m * o.n / pc, nil
}

// L computes L = ∂Cc/∂pc
func (o VanGen) L(pc, sl float64, wet bool) (float64, error) {
	if pc <= o.pcmin || pc >= o.pclim {
		return 0, nil
	}
	c := math.Pow(o.α*pc, o.n)
	mn := o.m * o.n
	fac := o.slmax - o.slmin
	return fac * c * math.Pow(c+1.0, -o.m-2.0) * mn * (c*mn - o.n + c + 1.0) / (pc * pc), nil
}

// J computes J = ∂Cc/∂sl
func (o VanGen) J(pc, sl float64, wet bool) (float64, error) {
	return 0, nil
}

// Derivs compute ∂Cc/∂pc and ∂²Cc/∂pc²
func (o VanGen) Derivs(pc, sl float64, wet bool) (L, Lx, J, Jx, Jy float64, err error) {
	if pc <= o.pcmin || pc >= o.pclim {
		return
	}
	c := math.Pow(o.α*pc, o.n)
	d := math.Pow(o.α*pc, o.n*2.0)
	mm := o.m * o.m
	nn := o.n * o.n
	mn := o.m * o.n
	ppp := pc * pc * pc
	fac := o.slmax - o.slmin
	L = fac * c * math.Pow(c+1.0, -o.m-2.0) * mn * (c*mn - o.n + c + 1.0) / (pc * pc)
	Lx = -fac * c * math.Pow(c+1.0, -o.m-3.0) * mn * (d*mm*nn - 3.0*c*o.m*nn - c*nn + nn + 3.0*d*mn + 3.0*c*mn - 3.0*c*o.n - 3.0*o.n + 2.0*d + 4.0*c + 2.0) / ppp
	return
}
