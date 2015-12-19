// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mreten

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/plt"
)

func Test_bc01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("bc01")

	mdl := new(BrooksCorey)
	prm := mdl.GetPrms(true)
	slmax := prm.Find("slmax")
	slmax.V = 0.95
	mdl.Init(prm)

	pc0 := -1.0
	sl0 := mdl.SlMax()
	pcf := 3.0
	nptsA := 11
	nptsB := 11

	if chk.Verbose {
		plt.Reset()
		Plot(mdl, pc0, sl0, pcf, nptsA, "'b.-'", "'r+-'", "bc")
	}

	tolCc := 1e-11
	tolD1a, tolD1b := 1e-10, 1e-17
	tolD2a, tolD2b := 1e-10, 1e-17
	Check(tst, mdl, pc0, sl0, pcf, nptsB, tolCc, tolD1a, tolD1b, tolD2a, tolD2b, chk.Verbose, []float64{0.2}, 1e-7, chk.Verbose)

	if chk.Verbose {
		PlotEnd(true)
	}
}
