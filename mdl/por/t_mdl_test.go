// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
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

func Test_mdl01(tst *testing.T) {

	//verbose()
	//doplot := true
	doplot := false
	chk.PrintTitle("mdl01")

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

	// constants
	H := 10.0
	grav := 10.0

	// fluids
	Liq := new(fld.Model)
	Liq.Init(Liq.GetPrms(true, false), H, grav)
	Gas := new(fld.Model)
	Gas.Init(Gas.GetPrms(true, true), H, grav)

	// porous model
	mdl := new(Model)
	err = mdl.Init(mdl.GetPrms(example), Cnd, Lrm, Liq, Gas, grav)
	if err != nil {
		tst.Errorf("Init failed: %v\n", err)
		return
	}
	//mdl.MEtrial = false
	mdl.ShowR = true

	// initial and final values
	pc0 := -5.0
	pcf := 20.0

	// plot lrm
	if doplot {
		npts := 41
		plt.Reset()
		lrm.Plot(mdl.Lrm, pc0, sl0, pcf, npts, "'b.-'", "", lrm_name)
	}

	// state A
	pcA := 5.0
	A, err := mdl.NewState(mdl.Liq.R0, mdl.Gas.R0, -pcA, 0)
	if err != nil {
		tst.Errorf("NewState failed: %v\n", err)
		return
	}

	// state B
	pcB := 10.0
	B, err := mdl.NewState(mdl.Liq.R0, mdl.Gas.R0, -pcB, 0)
	if err != nil {
		tst.Errorf("NewState failed: %v\n", err)
		return
	}

	// plot A and B points
	if doplot {
		plt.PlotOne(pcA, A.A_sl, "'gs', clip_on=0, label='A', ms=10")
		plt.PlotOne(pcB, B.A_sl, "'ks', clip_on=0, label='B'")
	}

	// incremental update
	test := 0
	var Δpl float64
	var n, iwet int
	switch test {
	case 1:
		Δpl = -5.0
		n, iwet = 23, 10
	case 2:
		Δpl = -20.0
		n, iwet = 10, 2
	default:
		Δpl = -1.0
		n, iwet = 41, 15
	}
	Pc := make([]float64, n)
	Sl := make([]float64, n)
	pl := -pcA
	Pc[0] = pcA
	Sl[0] = A.A_sl
	for i := 1; i < n; i++ {
		if i > iwet {
			Δpl = -Δpl
			iwet = n
		}
		pl += Δpl
		err = mdl.Update(A, Δpl, 0, pl, 0)
		if err != nil {
			tst.Errorf("test failed: %v\n", err)
			return
		}
		Pc[i] = -pl
		Sl[i] = A.A_sl
	}

	// show graph
	if doplot {
		plt.Plot(Pc, Sl, "'ro-', clip_on=0, label='update'")
		lrm.PlotEnd(true)
	}
}
