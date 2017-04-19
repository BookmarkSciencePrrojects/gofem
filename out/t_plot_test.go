// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"testing"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

func Test_plot01(tst *testing.T) {

	// test title
	//verbose()
	chk.PrintTitle("plot01")

	// start simulation
	main := fem.NewMain("data/p02.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// start post-processing
	Start("data/p02.sim", 0, 0)

	// define entities
	Define("A", N{1})
	Define("B", At{0, 1})
	//Define("left", Along{{0, 0}, {0, 10}})
	Define("left", AlongY{0}) // 0 => x_cte

	// load results
	LoadResults(nil)

	plA := GetRes("pl", "A", 0)

	Splot("t-pl", "liquid pressure")
	Plot("t", "pl", "B", &plt.A{C: "b", M: "."}, -1)
	Plot("t", plA, "A", &plt.A{C: "r", M: "."}, -1)

	Splot("pl-pl", "")
	Plot("pl", "pl", "A", &plt.A{C: "k", M: "o"}, -1)

	Splot("pl-y", "")
	Plot("pl", "y", "left", &plt.A{C: "b", M: "o"}, 0)
	Plot("pl", "y", "left", &plt.A{C: "g", M: "o"}, -1)

	Splot("y-pl", "")
	io.Pforan("T = %v\n", Times)
	last := len(Times) - 1
	Plot("y", "pl", "left", &plt.A{C: "b", M: "o", L: io.Sf("t=%g", Times[0])}, 0)
	Plot("y", "pl", "left", &plt.A{C: "m", M: "*", Lw: 2, L: io.Sf("t=%g", Times[last])}, -1)

	if chk.Verbose {
		Draw("", "", -1, -1, false, nil)
	}
}
