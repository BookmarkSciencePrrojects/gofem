// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
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
		sol.Init([]*fun.Prm{
			&fun.Prm{N: "r", V: 1.0},
			&fun.Prm{N: "E", V: 1e3},
			&fun.Prm{N: "nu", V: 0.3},
			&fun.Prm{N: "qnV", V: 0.0},
			&fun.Prm{N: "qnH", V: 10.0},
		})

		L := 4.0
		npts := 101

		plt.SetForEps(1.2, 455)
		sol.PlotStress(1, L, npts)
		plt.SaveD("/tmp/gofem", "ana_platehole01.eps")
	}
}

func Test_selfweight01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("selfweight01")

	var sol ConfinedSelfWeight
	sol.Init([]*fun.Prm{
		&fun.Prm{N: "E", V: 1e3},
		&fun.Prm{N: "nu", V: 0.25},
		&fun.Prm{N: "rho", V: 2.0},
		&fun.Prm{N: "g", V: 10.0},
		&fun.Prm{N: "h", V: 3.0},
		&fun.Prm{N: "w", V: 1.0},
	})

	x := []float64{0, 0, 0}
	σ := sol.Stress(1, x)
	io.Pforan("σ = %v\n", σ)
	chk.Scalar(tst, "σz @ z=0", 1e-17, σ[2], -60.0)

	x[2] = 1.5
	σ = sol.Stress(1, x)
	io.Pforan("σ = %v\n", σ)
	chk.Scalar(tst, "σz @ z=1.5", 1e-17, σ[2], -30.0)

	x[2] = 3.0
	u := sol.Displ(1, x)
	io.Pforan("u = %v\n", u)
}
