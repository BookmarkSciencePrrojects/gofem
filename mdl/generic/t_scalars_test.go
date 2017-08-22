// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package generic

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
)

func Test_scalars01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("scalars")

	mdl, err := New("scalars")
	if err != nil {
		tst.Errorf("New failed: %v\n", err)
		return
	}

	prms := []*fun.P{
		&fun.P{N: "A", V: 1},
		&fun.P{N: "B", V: 2},
		&fun.P{N: "C", V: 3},
		&fun.P{N: "D", V: 4},
		&fun.P{N: "E", V: 5},
		&fun.P{N: "F", V: 6},
		&fun.P{N: "G", V: 7},
		&fun.P{N: "H", V: 8},
		&fun.P{N: "I", V: 9},
		&fun.P{N: "J", V: 10},
		&fun.P{N: "K", V: 11},
		&fun.P{N: "L", V: 12},
		&fun.P{N: "M", V: 13},
		&fun.P{N: "N", V: 14},
		&fun.P{N: "O", V: 15},
		&fun.P{N: "P", V: 16},
		&fun.P{N: "Q", V: 17},
		&fun.P{N: "R", V: 18},
		&fun.P{N: "S", V: 19},
		&fun.P{N: "T", V: 20},
		&fun.P{N: "U", V: 21},
		&fun.P{N: "V", V: 22},
		&fun.P{N: "W", V: 23},
		&fun.P{N: "X", V: 24},
		&fun.P{N: "Y", V: 25},
		&fun.P{N: "Z", V: 26},
	}

	ndim := 3
	err = mdl.Init(ndim, prms)
	if err != nil {
		tst.Errorf("cannot initialise model: %v\n", err)
		return
	}

	m := mdl.(*Scalars)
	chk.Float64(tst, "A", 1e-15, m.A, 1)
	chk.Float64(tst, "B", 1e-15, m.B, 2)
	chk.Float64(tst, "C", 1e-15, m.C, 3)
	chk.Float64(tst, "D", 1e-15, m.D, 4)
	chk.Float64(tst, "E", 1e-15, m.E, 5)
	chk.Float64(tst, "F", 1e-15, m.F, 6)
	chk.Float64(tst, "G", 1e-15, m.G, 7)
	chk.Float64(tst, "H", 1e-15, m.H, 8)
	chk.Float64(tst, "I", 1e-15, m.I, 9)
	chk.Float64(tst, "J", 1e-15, m.J, 10)
	chk.Float64(tst, "K", 1e-15, m.K, 11)
	chk.Float64(tst, "L", 1e-15, m.L, 12)
	chk.Float64(tst, "M", 1e-15, m.M, 13)
	chk.Float64(tst, "N", 1e-15, m.N, 14)
	chk.Float64(tst, "O", 1e-15, m.O, 15)
	chk.Float64(tst, "P", 1e-15, m.P, 16)
	chk.Float64(tst, "Q", 1e-15, m.Q, 17)
	chk.Float64(tst, "R", 1e-15, m.R, 18)
	chk.Float64(tst, "S", 1e-15, m.S, 19)
	chk.Float64(tst, "T", 1e-15, m.T, 20)
	chk.Float64(tst, "U", 1e-15, m.U, 21)
	chk.Float64(tst, "V", 1e-15, m.V, 22)
	chk.Float64(tst, "W", 1e-15, m.W, 23)
	chk.Float64(tst, "X", 1e-15, m.X, 24)
	chk.Float64(tst, "Y", 1e-15, m.Y, 25)
	chk.Float64(tst, "Z", 1e-15, m.Z, 26)
}
