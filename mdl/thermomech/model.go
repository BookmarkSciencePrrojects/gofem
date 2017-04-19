// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package thermomech implements models thermo-mechanical problems
package thermomech

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
)

// Model defines the interface for thermomech models
type Model interface {
	Init(ndim int, prms fun.Prms) error // initialises model
}

// New thermomech model
func New(name string) (model Model, err error) {
	allocator, ok := allocators[name]
	if !ok {
		return nil, chk.Err("model %q is not available in 'thermomech' database", name)
	}
	return allocator(), nil
}

// allocators holds all available models
var allocators = map[string]func() Model{}
