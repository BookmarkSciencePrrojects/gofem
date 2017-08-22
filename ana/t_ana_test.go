// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package ana

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

func Test_platehole01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("platehole01")

	if chk.Verbose {

		var sol PlateHole
		sol.Init([]*fun.P{
			&fun.P{N: "r", V: 1.0},
			&fun.P{N: "E", V: 1e3},
			&fun.P{N: "nu", V: 0.3},
			&fun.P{N: "qnV", V: 0.0},
			&fun.P{N: "qnH", V: 10.0},
		})

		L := 4.0
		npts := 101

		plt.Reset(false, nil)
		sol.PlotStress(1, L, npts)
		plt.Save("/tmp/gofem", "ana_platehole01")
	}
}

func Test_selfweight01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("selfweight01")

	var sol ConfinedSelfWeight
	sol.Init([]*fun.P{
		&fun.P{N: "E", V: 1e3},
		&fun.P{N: "nu", V: 0.25},
		&fun.P{N: "rho", V: 2.0},
		&fun.P{N: "g", V: 10.0},
		&fun.P{N: "h", V: 3.0},
		&fun.P{N: "w", V: 1.0},
	})

	x := []float64{0, 0, 0}
	σ := sol.Stress(1, x)
	io.Pforan("σ = %v\n", σ)
	chk.Float64(tst, "σz @ z=0", 1e-17, σ[2], -60.0)

	x[2] = 1.5
	σ = sol.Stress(1, x)
	io.Pforan("σ = %v\n", σ)
	chk.Float64(tst, "σz @ z=1.5", 1e-17, σ[2], -30.0)

	x[2] = 3.0
	u := sol.Displ(1, x)
	io.Pforan("u = %v\n", u)
}
