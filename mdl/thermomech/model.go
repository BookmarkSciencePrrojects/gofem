// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package thermomech implements models thermo-mechanical problems
package thermomech

import "github.com/cpmech/gosl/fun"

// Model defines the interface for thermomech models
type Model interface {
	Init(ndim int, prms fun.Prms) error // initialises model
}

// allocators holds all available models
var allocators = map[string]func() Model{}
