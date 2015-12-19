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

func Test_mdl01(tst *testing.T) {

	//verbose()
	//doplot := true
	doplot := false
	chk.PrintTitle("mdl01")

	// info
	simfnk := "mdl01"
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
	//mdl.MEtrial = false
	mdl.ShowR = true

	// initial and final values
	pc0 := -5.0
	pcf := 20.0

	// plot lrm
	if doplot {
		npts := 41
		plt.Reset()
		mreten.Plot(mdl.Lrm, pc0, sl0, pcf, npts, "'b.-'", "", lrm_name)
	}

	// state A
	pcA := 5.0
	A, err := mdl.NewState(mdl.RhoL0, mdl.RhoG0, -pcA, 0)
	if err != nil {
		tst.Errorf("mporous.NewState failed: %v\n", err)
		return
	}

	// state B
	pcB := 10.0
	B, err := mdl.NewState(mdl.RhoL0, mdl.RhoG0, -pcB, 0)
	if err != nil {
		tst.Errorf("mporous.NewState failed: %v\n", err)
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
		mreten.PlotEnd(true)
	}
}
