// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"testing"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/tests"
	"github.com/cpmech/gosl/chk"
)

func Test_rjoint01(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("rjoint01. curved line in 3D")

	// initialisation
	main := fem.NewMain("data/rjoint01.sim", "", true, false, false, false, chk.Verbose, 0)

	// callback to check consistent tangent operators
	eid := 2   // rjoint element
	if false { // TODO: fix this
		tests.Rjoint(main, &tests.Kb{
			Tst: tst, Eid: eid, Tol: 1e-8, Verb: chk.Verbose,
			Ni: -1, Nj: -1, ItMin: 1, ItMax: -1, Tmin: 0.6, Tmax: 0.6,
		})
	}

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}
}
