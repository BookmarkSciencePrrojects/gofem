// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package gen implements a place holder for generic models
package gen

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
)

// Model defines a generic model type
type Model interface {
	Init(ndim int, prms fun.Prms) error // Init initialises this structure
}

// New generic model
func New(name string) (model Model, err error) {
	allocator, ok := allocators[name]
	if !ok {
		return nil, chk.Err("model %q is not available in 'gen' database", name)
	}
	return allocator(), nil
}

// allocators holds all available models
var allocators = map[string]func() Model{}
