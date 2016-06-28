// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package por

import (
	"testing"

	"github.com/cpmech/gofem/mdl/cnd"
	"github.com/cpmech/gofem/mdl/fld"
	"github.com/cpmech/gofem/mdl/lrm"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/plt"
)

func Test_plot01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("plot01")

	// conductivity model
	example := true
	Cnd := new(cnd.M1)
	err := Cnd.Init(Cnd.GetPrms(example))
	if err != nil {
		tst.Errorf("Init failed: %v\n", err)
		return
	}

	// liquid retention model
	lrm_name := "ref-m1"
	//lrm_name := "vg"
	var Lrm lrm.Model
	var prm fun.Prms
	sl0 := 0.95
	switch lrm_name {
	case "ref-m1":
		Lrm = new(lrm.RefM1)
		prm = Lrm.GetPrms(example)
		y0 := prm.Find("y0")
		y0.V = sl0
		nowet := prm.Find("nowet")
		nowet.V = 1
	case "vg":
		Lrm = new(lrm.VanGen)
		prm = Lrm.GetPrms(example)
		slmax := prm.Find("slmax")
		slmax.V = sl0
	}
	err = Lrm.Init(prm)
	if err != nil {
		tst.Errorf("nit failed: %v\n", err)
		return
	}

	// constants
	H := 10.0
	grav := 10.0

	// fluids
	Liq := new(fld.Model)
	Liq.Init(Liq.GetPrms(true), H, grav)
	Gas := new(fld.Model)
	Gas.Gas = true
	Gas.Init(Gas.GetPrms(true), H, grav)

	// porous model
	mdl := new(Model)
	err = mdl.Init(mdl.GetPrms(example), Cnd, Lrm, Liq, Gas, grav)
	if err != nil {
		tst.Errorf("Init failed: %v\n", err)
		return
	}
	mdl.ShowR = true

	// plot
	if chk.Verbose {
		np := 61
		plt.SetForEps(1.2, 350)
		PlotLrm(mdl, "/tmp/gofem", "fig_por_plot01.eps", 20, np, true, true, true, "", "", "", "", "")
	}
}
