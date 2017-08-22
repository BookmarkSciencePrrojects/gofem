// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package generic

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun/dbf"
)

func Test_scalars01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("scalars")

	mdl, err := New("scalars")
	if err != nil {
		tst.Errorf("New failed: %v\n", err)
		return
	}

	prms := []*dbf.P{
		&dbf.P{N: "A", V: 1},
		&dbf.P{N: "B", V: 2},
		&dbf.P{N: "C", V: 3},
		&dbf.P{N: "D", V: 4},
		&dbf.P{N: "E", V: 5},
		&dbf.P{N: "F", V: 6},
		&dbf.P{N: "G", V: 7},
		&dbf.P{N: "H", V: 8},
		&dbf.P{N: "I", V: 9},
		&dbf.P{N: "J", V: 10},
		&dbf.P{N: "K", V: 11},
		&dbf.P{N: "L", V: 12},
		&dbf.P{N: "M", V: 13},
		&dbf.P{N: "N", V: 14},
		&dbf.P{N: "O", V: 15},
		&dbf.P{N: "P", V: 16},
		&dbf.P{N: "Q", V: 17},
		&dbf.P{N: "R", V: 18},
		&dbf.P{N: "S", V: 19},
		&dbf.P{N: "T", V: 20},
		&dbf.P{N: "U", V: 21},
		&dbf.P{N: "V", V: 22},
		&dbf.P{N: "W", V: 23},
		&dbf.P{N: "X", V: 24},
		&dbf.P{N: "Y", V: 25},
		&dbf.P{N: "Z", V: 26},
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
