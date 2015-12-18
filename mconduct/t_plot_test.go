// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mconduct

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/plt"
)

func Test_plot01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("plot01")

	if !chk.Verbose {
		return
	}

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

	plt.SetForEps(1.2, 350)
	gas := false
	Plot(mdl, "gofem", "mconduct_plot01_liq.eps", 101, gas, true, true)

	plt.SetForEps(1.2, 350)
	gas = true
	Plot(mdl, "gofem", "mconduct_plot01_gas.eps", 101, gas, true, true)
}
