// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gosl/fun"
)

// Solver implements the actual solver (time loop)
type Solver interface {
	Run(tf float64, dtFunc, dtoFunc fun.Func, verbose bool, dbgKb DebugKb_t) (err error)
}

// allocators holds all available solvers
var allocators = make(map[string]func(doms []*Domain, sum *Summary, dc *ele.DynCoefs) Solver)
