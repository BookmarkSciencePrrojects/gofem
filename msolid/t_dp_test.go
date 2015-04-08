// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
)

func Test_dp01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("dp01")

	// allocate driver
	ndim, pstress := 2, false
	simfnk, modelname := "test", "dp"
	var drv Driver
	err := drv.Init(simfnk, modelname, ndim, pstress, []*fun.Prm{
		&fun.Prm{N: "K", V: 1.5},
		&fun.Prm{N: "G", V: 1},
		&fun.Prm{N: "M", V: 0},
		&fun.Prm{N: "Mb", V: 0},
		&fun.Prm{N: "qy0", V: 2},
		&fun.Prm{N: "H", V: 0.5},
	})
	drv.CheckD = true
	drv.VerD = false // verbose
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}

	// dp model
	dp := drv.model.(*DruckerPrager)

	// path
	p0 := 0.0
	Δp := 3.0
	Δq := dp.qy0 + dp.M*Δp
	ϵ := 1e-3
	DP := []float64{Δp + ϵ, 3, 2, 1, 0}
	DQ := []float64{Δq + ϵ, 4, 2, 1, 3}
	nincs := 1
	niout := 1
	noise := 0.0
	var pth Path
	err = pth.SetPQstrain(ndim, nincs, niout, dp.K, dp.G, p0, DP, DQ, noise)
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}

	// run
	err = drv.Run(&pth)
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}

	// plot
	//if true {
	if false {
		var plr Plotter
		plr.SetFig(false, false, 1, 400, "/tmp", "test_dp01")
		plr.SetModel(dp)
		plr.PreCor = drv.PreCor
		plr.Plot(PlotSet7, drv.Res, drv.Eps, true, true)
	}
}
