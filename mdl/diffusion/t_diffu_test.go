// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package diffusion

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun/dbf"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

func Test_m1(tst *testing.T) {

	//verbose()
	chk.PrintTitle("m1")

	mdl, err := New("m1")
	if err != nil {
		tst.Errorf("New failed: %v\n", err)
		return
	}

	prms := []*dbf.P{
		&dbf.P{N: "a0", V: 1.0},
		&dbf.P{N: "a1", V: 2.0},
		&dbf.P{N: "a2", V: 3.0},
		&dbf.P{N: "a3", V: 4.0},
		&dbf.P{N: "rho", V: 3.3},
		&dbf.P{N: "k", V: 0.1},
	}

	ndim := 3
	err = mdl.Init(ndim, prms)
	if err != nil {
		tst.Errorf("cannot initialise model: %v\n", err)
		return
	}

	m := mdl.(*M1)
	chk.Float64(tst, "a0", 1e-15, m.a0, 1.0)
	chk.Float64(tst, "a1", 1e-15, m.a1, 2.0)
	chk.Float64(tst, "a2", 1e-15, m.a2, 3.0)
	chk.Float64(tst, "a3", 1e-15, m.a3, 4.0)
	chk.Float64(tst, "rho", 1e-15, m.Rho, 3.3)
	chk.Deep2(tst, "kcte", 1e-15, m.Kcte, [][]float64{
		{0.1, 0, 0},
		{0, 0.1, 0},
		{0, 0, 0.1},
	})

	u := 0.5
	kten := la.MatAlloc(ndim, ndim)
	m.Kten(kten, u)
	kval := 1.0 + 2.0*u + 3.0*u*u + 4.0*u*u*u
	chk.Float64(tst, "kval", 1e-15, m.Kval(u), kval)
	chk.Deep2(tst, "kcte", 1e-15, kten, [][]float64{
		{0.1 * kval, 0, 0},
		{0, 0.1 * kval, 0},
		{0, 0, 0.1 * kval},
	})

	U := utl.LinSpace(0, 2.0, 5)
	for _, uval := range U {
		dana := m.DkDu(uval)
		chk.DerivScaSca(tst, "DkDu", 1e-9, dana, uval, 1e-3, chk.Verbose, func(x float64) (float64, error) {
			return m.Kval(x), nil
		})
	}
}
