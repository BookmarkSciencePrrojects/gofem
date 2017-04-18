// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"math"
	"sort"
	"testing"

	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/ele/solid"
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/tests"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func Test_beam01a(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("beam01a. check DOFs")

	// fem
	main := fem.NewMain("data/beam01.sim", "", true, false, false, false, chk.Verbose, 0)

	// set stage
	err := main.SetStage(0)
	if err != nil {
		tst.Errorf("SetStage failed:\n%v", err)
		return
	}

	// initialise solution vectors
	err = main.ZeroStage(0, true)
	if err != nil {
		tst.Errorf("ZeroStage failed:\n%v", err)
		return
	}

	// domain
	dom := main.Domains[0]

	// nodes and elements
	chk.IntAssert(len(dom.Nodes), 2)
	chk.IntAssert(len(dom.Elems), 1)

	// check dofs
	for _, nod := range dom.Nodes {
		chk.IntAssert(len(nod.Dofs), 3)
	}

	// check equations
	nids, eqs := tests.GetNidsEqs(dom)
	chk.Ints(tst, "nids", nids, []int{0, 1})
	chk.Ints(tst, "eqs", eqs, []int{0, 1, 2, 3, 4, 5})

	// check solution arrays
	ny := 6
	nλ := 3
	nyb := ny + nλ
	chk.IntAssert(len(dom.Sol.Y), ny)
	chk.IntAssert(len(dom.Sol.Dydt), 0)
	chk.IntAssert(len(dom.Sol.D2ydt2), 0)
	chk.IntAssert(len(dom.Sol.Psi), 0)
	chk.IntAssert(len(dom.Sol.Zet), 0)
	chk.IntAssert(len(dom.Sol.Chi), 0)
	chk.IntAssert(len(dom.Sol.L), nλ)
	chk.IntAssert(len(dom.Sol.ΔY), ny)

	// check linear solver arrays
	chk.IntAssert(len(dom.Fb), nyb)
	chk.IntAssert(len(dom.Wb), nyb)

	// check umap
	e := dom.Elems[0].(*solid.Beam)
	io.Pforan("e = %v\n", e.Umap)
	chk.Ints(tst, "umap", e.Umap, []int{0, 1, 2, 3, 4, 5})

	// constraints
	chk.IntAssert(len(dom.EssenBcs.Bcs), nλ)
	var ct_ux_eqs []int // constrained ux equations [sorted]
	var ct_uy_eqs []int // constrained uy equations [sorted]
	for _, c := range dom.EssenBcs.Bcs {
		chk.IntAssert(len(c.Eqs), 1)
		eq := c.Eqs[0]
		io.Pforan("key=%v eq=%v\n", c.Key, eq)
		switch c.Key {
		case "ux":
			ct_ux_eqs = append(ct_ux_eqs, eq)
		case "uy":
			ct_uy_eqs = append(ct_uy_eqs, eq)
		default:
			tst.Errorf("key %s is incorrect", c.Key)
			return
		}
	}
	sort.Ints(ct_ux_eqs)
	sort.Ints(ct_uy_eqs)
	chk.Ints(tst, "constrained ux equations", ct_ux_eqs, []int{0})
	chk.Ints(tst, "constrained uy equations", ct_uy_eqs, []int{1, 4})
}

func Test_beam01b(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("beam01b. simply supported")

	// start simulation
	main := fem.NewMain("data/beam01.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	dom := main.Domains[0]
	e := dom.Elems[0].(*solid.Beam)
	M22 := e.CalcMoment2d(dom.Sol, 0.5, 1)
	qn, L := 15.0, 1.0
	Mcentre := qn * L * L / 8.0
	io.Pforan("M22 = %v (%v)\n", M22, Mcentre)
	chk.Scalar(tst, "M22 @ centre", 1e-17, M22[0], Mcentre)

	// check moment using OutIpsData
	idx_centre := 5 // considering 11 stations
	res := ele.NewIpsMap()
	e.OutIpVals(res, dom.Sol)
	io.Pfcyan("M22 @ centre (OutIpsData) = %v\n", res.Get("M22", idx_centre))
	chk.Scalar(tst, "M22 @ centre (OutIpsData)", 1e-17, res.Get("M22", idx_centre), Mcentre)
}

func Test_beam02(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("beam02. cantilever")

	// start simulation
	main := fem.NewMain("data/beam02.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	dom := main.Domains[0]
	ele := dom.Elems[0].(*solid.Beam)
	M22 := ele.CalcMoment2d(dom.Sol, 0, 1)
	qn, L := 15.0, 1.0
	Mleft := -qn * L * L / 2.0
	io.Pforan("M22 = %v (%v)\n", M22, Mleft)
	chk.Scalar(tst, "M @ left", 1e-15, M22[0], Mleft)
}

func Test_beam03(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("beam03. small frame")

	// start simulation
	main := fem.NewMain("data/beam03.sim", "", true, true, false, false, chk.Verbose, 0)

	// set stage and run
	set_and_run := func(stgidx int) {
		err := main.SetStage(stgidx)
		if err != nil {
			tst.Errorf("SetStage failed:\n%v", err)
			return
		}
		err = main.SolveOneStage(stgidx, true)
		if err != nil {
			tst.Error("SolveOneStage failed:\n%v", err)
			return
		}
	}
	set_and_run(0)

	// domain
	dom := main.Domains[0]

	// message
	io.Pf("\nebc = %v\n", dom.EssenBcs.List(0))
	io.Pf("fbc = %v\n\n", dom.PtNatBcs.List(0))
	for _, nod := range dom.Nodes {
		io.Pf("node # %2d ux=%23.15e uy=%23.15e\n", nod.Vert.Id, dom.Sol.Y[nod.GetEq("ux")], dom.Sol.Y[nod.GetEq("uy")])
	}
	io.Pf("\n")

	// define function to check bending moment
	check_M := func(beamId int, ξ, Mref, tol float64) {
		ele := dom.Cid2elem[beamId].(*solid.Beam)
		M22 := ele.CalcMoment2d(dom.Sol, ξ, 1)
		chk.Scalar(tst, io.Sf("Beam %d: M22(ξ=%g) = %.6f", ele.Id(), ξ, M22[0]), tol, M22[0], Mref)
	}

	// check
	check_M(0, 0, 0, 1e-13)
	check_M(0, 1, 34, 1e-13)
	check_M(1, 0, -20.4, 1e-13)
	check_M(1, 1, 0, 1e-13)
	check_M(2, 0, 0, 1e-13)
	check_M(2, 1, -54.4, 1e-13)

	//nstations, withtext, numfmt, tol, coef, sf, onlyLin := 11, true, "", 1e-10, 0.2, 0.0, false
	if chk.Verbose {
		//plt.SetForPng(1, 600, 150)
		//dom.Msh.Draw2d(onlyLin, true, nil, 3)
		//fem.PlotAllBendingMom2d(dom, nstations, withtext, numfmt, tol, coef, sf)
		//plt.SaveD("/tmp/gofem", "test_beam03_prob1.png")
	}

	// problem # 2
	set_and_run(1)
	check_M(0, 0, 0, 1e-13)
	check_M(0, 1, 34, 1e-13)
	check_M(1, 0, -20.4, 1e-13)
	check_M(1, 1, 0, 1e-13)
	check_M(2, 0, 0, 1e-13)
	check_M(2, 1, 0, 1e-13)

	if chk.Verbose {
		//plt.SetForPng(1, 600, 150)
		//dom.Msh.Draw2d(onlyLin, true, nil, 3)
		//fem.PlotAllBendingMom2d(dom, nstations, withtext, numfmt, tol, coef, sf)
		//plt.SaveD("/tmp/gofem", "test_beam03_prob2.png")
	}

	// problem # 3
	set_and_run(2)
	check_M(0, 0, 0, 1e-13)
	check_M(0, 1, 20, 1e-12)
	check_M(1, 0, 20, 1e-12)
	check_M(1, 1, 0, 1e-13)
	check_M(2, 0, 0, 1e-13)
	check_M(2, 1, 0, 1e-13)

	if chk.Verbose {
		//plt.SetForPng(1, 600, 150)
		//dom.Msh.Draw2d(onlyLin, true, nil, 3)
		//fem.PlotAllBendingMom2d(dom, nstations, withtext, numfmt, tol, coef, sf)
		//plt.SaveD("/tmp/gofem", "test_beam03_prob3.png")
	}

	// problem # 4
	set_and_run(3)
	e0M := func(x float64) float64 { return 5.0752*x - 0.9984*x*x/3.0 - (16.0-x)*0.0624*x*x/6.0 }
	e1M := func(x float64) float64 { return 1.7472*(10.0-x) - 0.6*0.0624*math.Pow(10.0-x, 3.0)/6.0 }
	check_M(0, 0, e0M(0), 1e-13)
	check_M(0, 0.5, e0M(5), 1e-13)
	check_M(0, 1, e0M(10), 1e-12)
	check_M(1, 0, e1M(0), 1e-12)
	check_M(1, 0.5, e1M(5), 1e-12)
	check_M(1, 1, e1M(10), 1e-13)
	check_M(2, 0, 0, 1e-13)
	check_M(2, 1, 0, 1e-13)

	if chk.Verbose {
		//plt.SetForPng(1, 600, 150)
		//dom.Msh.Draw2d(onlyLin, true, nil, 3)
		//fem.PlotAllBendingMom2d(dom, nstations, withtext, numfmt, tol, coef, sf)
		//plt.SaveD("/tmp/gofem", "test_beam03_prob4.png")
	}
}

func Test_beam04(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("beam04. 3D beam (bh414)")

	// fem
	main := fem.NewMain("data/beam04.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 1e-10
	tolu := 1e-16
	tols := 1e-12
	tests.CompareResults(tst, "data/beam04.sim", "cmp/bh414.cmp", "", tolK, tolu, tols, skipK, chk.Verbose, nil)
}

func Test_beam05(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("beam04. 3D frame")

	// fem
	main := fem.NewMain("data/frame01.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// displacement @ top
	dom := main.Domains[0]
	nod := dom.Vid2node[7]
	eqx := nod.GetEq("ux")
	eqy := nod.GetEq("uy")
	eqz := nod.GetEq("uz")
	u := []float64{dom.Sol.Y[eqx], dom.Sol.Y[eqy], dom.Sol.Y[eqz]}
	io.Pforan("u @ sel node = %v\n", u)

	// TODO: add tests here
}
