// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package mconduct implements models for liquid and gas conductivity in porous media
package mconduct

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
)

// Model defines liquid-gas conductivity models
type Model interface {
	Init(prms fun.Prms) error      // Init initialises this structure
	GetPrms(example bool) fun.Prms // gets (an example) of parameters
	Klr(sl float64) float64        // Klr returns klr
	Kgr(sg float64) float64        // Kgr returns kgr
	DklrDsl(sl float64) float64    // DklrDsl returns ∂klr/∂sl
	DkgrDsg(sg float64) float64    // DkgrDsl returns ∂kgr/∂sl
}

// New conductivity model
func New(name string) (model Model, err error) {
	allocator, ok := allocators[name]
	if !ok {
		return nil, chk.Err("model %q is not available in mconduct database", name)
	}
	return allocator(), nil
}

// allocators holds all available models
var allocators = map[string]func() Model{}
