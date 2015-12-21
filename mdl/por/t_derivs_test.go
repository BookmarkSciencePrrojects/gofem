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
	Cnd := new(cnd.M1)
	err := Cnd.Init(Cnd.GetPrms(example))
	if err != nil {
		tst.Errorf("Init failed: %v\n", err)
		return
	}

	// liquid retention model
	lrm_name := "ref-m1"
	Lrm := new(lrm.RefM1)
	err = Lrm.Init(Lrm.GetPrms(example))
	if err != nil {
		tst.Errorf("Init failed: %v\n", err)
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
		plt.Reset()
		lrm.Plot(mdl.Lrm, pc0, 1.0, pcf, npts, "'b.-'", "'r+-'", lrm_name)
		n := len(drv.Res)
		Sl := make([]float64, n)
		for i, s := range drv.Res {
			Sl[i] = s.A_sl
		}
		plt.Plot(Pc, Sl, "'ko--', clip_on=0")
		lrm.PlotEnd(true)
	}
}
