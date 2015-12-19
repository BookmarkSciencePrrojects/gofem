// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mporous

import (
	"testing"

	"github.com/cpmech/gofem/mconduct"
	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/plt"
)

func Test_plot01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("plot01")

	if !chk.Verbose {
		return
	}

	// info
	simfnk := "derivs01"
	matname := "mat1"
	getnew := false
	example := true

	// conductivity model
	cnd := mconduct.GetModel(simfnk, matname, "m1", getnew)
	err := cnd.Init(cnd.GetPrms(example))
	if err != nil {
		tst.Errorf("mconduct.Init failed: %v\n", err)
		return
	}

	// liquid retention model
	lrm_name := "ref-m1"
	//lrm_name := "vg"
	lrm := mreten.GetModel(simfnk, matname, lrm_name, getnew)
	prm := lrm.GetPrms(example)
	sl0 := 0.95
	switch lrm_name {
	case "ref-m1":
		y0 := prm.Find("y0")
		y0.V = sl0
	case "vg":
		slmax := prm.Find("slmax")
		slmax.V = sl0
	}
	err = lrm.Init(prm)
	if err != nil {
		tst.Errorf("mreten.Init failed: %v\n", err)
		return
	}

	// porous model
	mdl := GetModel(simfnk, matname, getnew)
	err = mdl.Init(mdl.GetPrms(example), cnd, lrm)
	if err != nil {
		tst.Errorf("mporous.Init failed: %v\n", err)
		return
	}

	// plot
	plt.SetForEps(1.2, 350)
	PlotSimple(mdl, "gofem", "mporous_plot01.eps", 20, 101, true, true, true)
}
