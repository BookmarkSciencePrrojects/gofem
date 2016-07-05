// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package retention

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/plt"
)

func Test_refm1a(tst *testing.T) {

	//verbose()
	chk.PrintTitle("refm1a")

	mdl := new(RefM1)
	prm := mdl.GetPrms(true)
	y0 := prm.Find("y0")
	y0.V = 0.95
	nowet := prm.Find("nowet")
	nowet.V = 1
	mdl.Init(prm)

	pc0 := -5.0
	sl0 := mdl.SlMax()
	pcf := 20.0
	nptsA := 21
	nptsB := 21

	if chk.Verbose {
		plt.Reset()
		Plot(mdl, pc0, sl0, pcf, nptsA, false, "'b.-',markevery=10", "'r+-'", "ref-m1_drying")
	}

	tolCc := 1e-17
	tolD1a, tolD1b := 1e-10, 1e-10
	tolD2a, tolD2b := 1e-10, 1e-8
	Check(tst, mdl, pc0, sl0, pcf, nptsB, tolCc, tolD1a, tolD1b, tolD2a, tolD2b, chk.Verbose, []float64{0}, 1e-7, false)

	slf, err := Update(mdl, pc0, sl0, pcf-pc0)
	if err != nil {
		tst.Errorf("update failed: %v\n", err)
		return
	}

	if chk.Verbose {
		Plot(mdl, pcf, slf, pc0, nptsA, false, "'b*-',markevery=10", "'r+:'", "ref-m1_wetting")
	}

	tolD1b = 1e-4
	tolD2a, tolD2b = 1e-10, 1e-8
	Check(tst, mdl, pcf, slf, pc0, nptsB, tolCc, tolD1a, tolD1b, tolD2a, tolD2b, chk.Verbose, []float64{0}, 1e-7, false)

	if chk.Verbose {
		PlotEnd(false)
		plt.SaveD("/tmp/gofem", "fig_refm1a.eps")
	}
}
