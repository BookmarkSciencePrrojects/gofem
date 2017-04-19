// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package conduct

import (
	"testing"

	"github.com/cpmech/gosl/chk"
)

func Test_plot01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("plot01")

	mdl := new(M1)
	prm := mdl.GetPrms(true)

	betg := prm.Find("betg")
	betg.V = 20.0

	lam1g := prm.Find("lam1g")
	lam1g.V = 0.1

	alpg := prm.Find("alpg")
	alpg.V = 0.0

	err := mdl.Init(prm)
	if err != nil {
		tst.Errorf("Init failed:\n%v")
		return
	}

	if chk.Verbose {
		gas := false
		Plot(mdl, "/tmp/gofem", "cnd_plot01_liq.eps", 101, gas, true, true)

		gas = true
		Plot(mdl, "/tmp/gofem", "cnd_plot01_gas.eps", 101, gas, true, true)
	}
}
