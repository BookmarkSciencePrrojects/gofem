// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package por

import (
	"testing"

	"github.com/cpmech/gofem/mdl/cnd"
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

	// porous model
	mdl := new(Model)
	err = mdl.Init(mdl.GetPrms(example), Cnd, Lrm)
	if err != nil {
		tst.Errorf("Init failed: %v\n", err)
		return
	}

	// plot
	if chk.Verbose {
		plt.SetForEps(1.2, 350)
		PlotSimple(mdl, "/tmp/gofem", "fig_por_plot01.eps", 20, 101, true, true, true)
	}
}
