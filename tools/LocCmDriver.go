// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"encoding/json"
	"path"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mdl/solid"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/tsr"
)

type Input struct {
	Dir     string
	SimFn   string
	MatName string
	PathFn  string
	PlotSet []string
	FigEps  bool
	FigProp float64
	FigWid  float64

	// derived
	inpfn string
}

func (o *Input) PostProcess() {
	if len(o.PlotSet) == 0 {
		o.PlotSet = solid.PlotSet6
	}
	if o.FigProp < 0.1 {
		o.FigProp = 1.0
	}
	if o.FigWid < 10 {
		o.FigWid = 400
	}
}

func (o Input) String() (l string) {
	l = io.ArgsTable("INPUT ARGUMENTS",
		"input filename", "inpfn", o.inpfn,
		"directory with .sim and .pat files", "Dir", o.Dir,
		"simulation filename", "SimFn", o.SimFn,
		"material name", "MatName", o.MatName,
		"path filename", "PathFn", o.PathFn,
		"plot set", "PlotSet", io.Sf("%v", o.PlotSet),
		"fig: generate .eps instead of .png", "FigEps", o.FigEps,
		"fig: proportion of figure", "FigProp", o.FigProp,
		"fig: width  of figure", "FigWid", o.FigWid,
	)
	return
}

func main() {

	// catch errors
	defer func() {
		if err := recover(); err != nil {
			io.PfRed("ERROR: %v\n", err)
		}
	}()

	// input data file
	var in Input
	in.inpfn, _ = io.ArgToFilename(0, "data/loccmdrv1", ".inp", true)

	// read and parse input data
	b, err := io.ReadFile(in.inpfn)
	if err != nil {
		io.PfRed("cannot read %s\n", in.inpfn)
		return
	}
	err = json.Unmarshal(b, &in)
	if err != nil {
		io.PfRed("cannot parse %s\n", in.inpfn)
		return
	}
	in.PostProcess()

	// print input table
	io.Pf("%v\n", in)

	// load simulation
	sim := inp.ReadSim(in.Dir+"/"+in.SimFn, "", false, 0)
	if sim == nil {
		io.PfRed("cannot load simulation\n")
		return
	}

	// get material data
	mat := sim.MatModels.Get(in.MatName)
	if mat == nil {
		io.PfRed("cannot get material\n")
		return
	}
	mdl := mat.Sld
	if mdl == nil {
		io.PfRed("cannot get solid model\n")
		return
	}

	// constants
	ndim := 3

	// load path
	var pth solid.Path
	err = pth.ReadJson(ndim, path.Join(in.Dir, in.PathFn))
	if err != nil {
		io.PfRed("cannot read path file %v\n", err)
		return
	}
	//io.PfYel("pth = %v\n", pth)

	// driver
	var drv solid.Driver
	drv.InitWithModel(ndim, mdl)

	// run
	err = drv.Run(&pth)
	if err != nil {
		io.Pfred("driver: Run failed: %v\n", err)
	}

	// plot
	//if false {
	if true {
		var plr solid.Plotter
		plr.SetFig(false, in.FigEps, in.FigProp, in.FigWid, "/tmp", "cmd_"+in.SimFn)
		var epm solid.EPmodel
		if m, ok := mdl.(solid.EPmodel); ok {
			plr.SetModel(m)
			epm = m
		}
		if epm != nil {
			//plr.Phi = epm.Get_phi()
			b := epm.Get_bsmp()
			epm.Set_bsmp(0)
			plr.YsClr0 = "magenta"
			plr.Plot(in.PlotSet, drv.Res, nil, true, false)
			epm.Set_bsmp(b)
		}
		plr.YsClr0 = "green"
		plr.Plot(in.PlotSet, drv.Res, drv.Eps, false, true)
	}

	// plot ys
	if false {
		//if true {
		plt.Reset()
		m := mdl.(*solid.SmpInvs)
		φ := m.Get_phi()
		σcCte := 10.0
		M := tsr.Phi2M(φ, "oct")
		rmin, rmax := 0.0, 1.28*M*σcCte
		nr, nα := 31, 81
		//nr,   nα   := 31, 1001
		npolarc := true
		simplec := false
		only0 := false
		grads := false
		showpts := false
		ferr := 10.0
		tsr.PlotOct("fig_isofun02.png", σcCte, rmin, rmax, nr, nα, φ, m.Isof.Fa, m.Isof.Ga,
			npolarc, simplec, only0, grads, showpts, true, true, ferr)
	}
}
