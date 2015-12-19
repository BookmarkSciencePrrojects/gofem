// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mporous

import (
	"testing"

	"github.com/cpmech/gofem/mconduct"
	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/plt"
)

func Test_plot01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("plot01")

	if !chk.Verbose {
		return
	}

	// conductivity model
	example := true
	cnd := new(mconduct.M1)
	err := cnd.Init(cnd.GetPrms(example))
	if err != nil {
		tst.Errorf("mconduct.Init failed: %v\n", err)
		return
	}

	// liquid retention model
	lrm_name := "ref-m1"
	//lrm_name := "vg"
	var lrm mreten.Model
	var prm fun.Prms
	sl0 := 0.95
	switch lrm_name {
	case "ref-m1":
		lrm = new(mreten.RefM1)
		prm = lrm.GetPrms(example)
		y0 := prm.Find("y0")
		y0.V = sl0
	case "vg":
		lrm = new(mreten.VanGen)
		prm = lrm.GetPrms(example)
		slmax := prm.Find("slmax")
		slmax.V = sl0
	}
	err = lrm.Init(prm)
	if err != nil {
		tst.Errorf("mreten.Init failed: %v\n", err)
		return
	}

	// porous model
	mdl := new(Model)
	err = mdl.Init(mdl.GetPrms(example), cnd, lrm)
	if err != nil {
		tst.Errorf("mporous.Init failed: %v\n", err)
		return
	}

	// plot
	plt.SetForEps(1.2, 350)
	PlotSimple(mdl, "gofem", "mporous_plot01.eps", 20, 101, true, true, true)
}
