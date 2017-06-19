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
	chk.Scalar(tst, "A", 1e-15, m.A, 1)
	chk.Scalar(tst, "B", 1e-15, m.B, 2)
	chk.Scalar(tst, "C", 1e-15, m.C, 3)
	chk.Scalar(tst, "D", 1e-15, m.D, 4)
	chk.Scalar(tst, "E", 1e-15, m.E, 5)
	chk.Scalar(tst, "F", 1e-15, m.F, 6)
	chk.Scalar(tst, "G", 1e-15, m.G, 7)
	chk.Scalar(tst, "H", 1e-15, m.H, 8)
	chk.Scalar(tst, "I", 1e-15, m.I, 9)
	chk.Scalar(tst, "J", 1e-15, m.J, 10)
	chk.Scalar(tst, "K", 1e-15, m.K, 11)
	chk.Scalar(tst, "L", 1e-15, m.L, 12)
	chk.Scalar(tst, "M", 1e-15, m.M, 13)
	chk.Scalar(tst, "N", 1e-15, m.N, 14)
	chk.Scalar(tst, "O", 1e-15, m.O, 15)
	chk.Scalar(tst, "P", 1e-15, m.P, 16)
	chk.Scalar(tst, "Q", 1e-15, m.Q, 17)
	chk.Scalar(tst, "R", 1e-15, m.R, 18)
	chk.Scalar(tst, "S", 1e-15, m.S, 19)
	chk.Scalar(tst, "T", 1e-15, m.T, 20)
	chk.Scalar(tst, "U", 1e-15, m.U, 21)
	chk.Scalar(tst, "V", 1e-15, m.V, 22)
	chk.Scalar(tst, "W", 1e-15, m.W, 23)
	chk.Scalar(tst, "X", 1e-15, m.X, 24)
	chk.Scalar(tst, "Y", 1e-15, m.Y, 25)
	chk.Scalar(tst, "Z", 1e-15, m.Z, 26)
}
