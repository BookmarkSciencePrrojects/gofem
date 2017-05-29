// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package porous

import (
	"testing"

	"github.com/cpmech/gofem/mdl/conduct"
	"github.com/cpmech/gofem/mdl/fluid"
	"github.com/cpmech/gofem/mdl/retention"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

func Test_derivs01(tst *testing.T) {

	//verbose()
	//doplot := true
	doplot := false
	chk.PrintTitle("derivs01")

	// conductivity model
	example := true
	Cnd := new(conduct.M1)
	err := Cnd.Init(Cnd.GetPrms(example))
	if err != nil {
		tst.Errorf("Init failed: %v\n", err)
		return
	}

	// liquid retention model
	lrm_name := "ref-m1"
	Lrm := new(retention.RefM1)
	err = Lrm.Init(Lrm.GetPrms(example))
	if err != nil {
		tst.Errorf("Init failed: %v\n", err)
		return
	}

	// constants
	H := 10.0
	grav := 10.0

	// fluids
	Liq := new(fluid.Model)
	Liq.Init(Liq.GetPrms(true), H, grav)
	Gas := new(fluid.Model)
	Gas.Gas = true
	Gas.Init(Gas.GetPrms(true), H, grav)

	// porous model
	mdl := new(Model)
	err = mdl.Init(mdl.GetPrms(example), Cnd, Lrm, Liq, Gas, grav)
	if err != nil {
		tst.Errorf("Init failed: %v\n", err)
		return
	}
	//mdl.Ncns = true
	//mdl.Ncns2 = true

	// path
	pc0 := 0.0
	pcf := 20.0
	np := 5
	//P := []float64{10, 5, 20, 0}
	P := []float64{5}
	Pc := GetPathCycle(pc0, P, np)
	io.Pforan("Pc = %v\n", Pc)

	// driver
	var drv Driver
	drv.TstD = tst
	err = drv.Init(mdl)
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}
	err = drv.Run(Pc)
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}

	// plot
	if doplot {
		npts := 41
		plt.Reset(false, nil)
		retention.Plot(mdl.Lrm, pc0, 1.0, pcf, npts, false, "'b.-'", "'r+-'", lrm_name)
		n := len(drv.Res)
		Sl := make([]float64, n)
		for i, s := range drv.Res {
			Sl[i] = s.A_sl
		}
		//plt.Plot(Pc, Sl, "'ko--', clip_on=0")
	}
}
