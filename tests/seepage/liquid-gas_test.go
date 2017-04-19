// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"sort"
	"testing"

	"github.com/cpmech/gofem/ele/seepage"
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/tests"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

func Test_pp01a(tst *testing.T) {

	/* this tests simulates seepage flow along a column
	 * by reducing the initial hydrostatic pressure at
	 * at the bottom of the column -- liquid-gas version
	 *
	 *      Nodes / Tags                    Equations
	 *
	 *     8 o----o----o 9 (-5)      44 45 o----o----o 42 43
	 *       |   14    |                   |  48 49  |
	 *       |         |                   |         |
	 *    21 o    o    o 22 (-6)     50 51 o    o    o 46 47
	 *       |   26    |                   |  52 53  |
	 *       |         |                   |         |
	 *     6 o----o----o 7 (-4)      32 33 o----o----o 30 31
	 *       |   13    |                   |  36 37  |
	 *       |         |                   |         |
	 *    19 o    o    o 20 (-6)     38 39 o    o    o 34 35
	 *       |   25    |                   |  40 41  |
	 *       |         |                   |         |
	 *     4 o----o----o 5 (-3)      20 21 o----o----o 18 19
	 *       |   12    |                   |  24 25  |
	 *       |         |                   |         |
	 *    17 o    o    o 18 (-6)     26 27 o    o    o 22 23
	 *       |   24    |                   |  28 29  |
	 *       |         |                   |         |
	 *     2 o----o----o 3 (-2)       6  7 o----o----o  4  5
	 *       |   11    |                   |  12 13  |
	 *       |         |                   |         |
	 *    15 o    o    o 16 (-6)     14 15 o    o    o 10 11
	 *       |   23    |                   |  16 17  |
	 *       |         |                   |         |
	 *     0 o----o----o 1 (-1)       0  1 o----o----o  2  3
	 *           10                            8 9
	 */

	//tests.Verbose()
	chk.PrintTitle("pp01a. Liquid-Gas. Check DOFs and BCs")

	// start simulation
	main := fem.NewMain("data/pp01.sim", "", true, false, false, false, chk.Verbose, 0)

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
	chk.IntAssert(len(dom.Nodes), 27)
	chk.IntAssert(len(dom.Elems), 4)

	// check dofs
	for _, nod := range dom.Nodes {
		chk.IntAssert(len(nod.Dofs), 2)
		chk.StrAssert(nod.Dofs[0].Key, "pl")
		chk.StrAssert(nod.Dofs[1].Key, "pg")
	}

	// check equations
	nids, eqs := tests.GetNidsEqs(dom)
	chk.Ints(tst, "nids", nids, []int{
		0, 1, 3, 2, 10, 16, 11, 15, 23,
		5, 4, 18, 12, 17, 24,
		7, 6, 20, 13, 19, 25,
		9, 8, 22, 14, 21, 26,
	})
	chk.Ints(tst, "eqs", eqs, utl.IntRange(54))

	// check pmap
	Plmaps := [][]int{
		{0, 2, 4, 6, 8, 10, 12, 14, 16},
		{6, 4, 18, 20, 12, 22, 24, 26, 28},
		{20, 18, 30, 32, 24, 34, 36, 38, 40},
		{32, 30, 42, 44, 36, 46, 48, 50, 52},
	}
	Pgmaps := [][]int{
		{1, 3, 5, 7, 9, 11, 13, 15, 17},
		{7, 5, 19, 21, 13, 23, 25, 27, 29},
		{21, 19, 31, 33, 25, 35, 37, 39, 41},
		{33, 31, 43, 45, 37, 47, 49, 51, 53},
	}
	for i, ele := range dom.Elems {
		e := ele.(*seepage.LiquidGas)
		io.Pforan("e%d.Plmap = %v\n", e.Id(), e.Plmap)
		io.Pforan("e%d.Pgmap = %v\n", e.Id(), e.Pgmap)
		chk.Ints(tst, "Plmap", e.Plmap, Plmaps[i])
		chk.Ints(tst, "Pgmap", e.Pgmap, Pgmaps[i])
	}

	// constraints
	chk.IntAssert(len(dom.EssenBcs.Bcs), 9)
	var ct_pl_eqs []int // equations with pl prescribed [sorted]
	var ct_pg_eqs []int // equations with pg prescribed [sorted]
	for _, c := range dom.EssenBcs.Bcs {
		chk.IntAssert(len(c.Eqs), 1)
		eq := c.Eqs[0]
		io.Pforan("key=%v eq=%v\n", c.Key, eq)
		switch c.Key {
		case "pl":
			ct_pl_eqs = append(ct_pl_eqs, eq)
		case "pg":
			ct_pg_eqs = append(ct_pg_eqs, eq)
		default:
			tst.Errorf("key %s is incorrect", c.Key)
		}
	}
	sort.Ints(ct_pl_eqs)
	sort.Ints(ct_pg_eqs)
	chk.Ints(tst, "equations with pl prescribed", ct_pl_eqs, []int{0, 2, 8})
	chk.Ints(tst, "equations with pg prescribed", ct_pg_eqs, []int{1, 3, 9, 43, 45, 49})

	// fluid models
	if dom.Sim.LiqMdl == nil || dom.Sim.GasMdl == nil {
		tst.Errorf("both liquid and gas materials must be specified in 'iniporous'\n")
		return
	}

	// initial values @ nodes
	io.Pforan("initial values @ nodes\n")
	for _, nod := range dom.Nodes {
		z := nod.Vert.C[1]
		eql := nod.Dofs[0].Eq
		eqg := nod.Dofs[1].Eq
		pl := dom.Sol.Y[eql]
		pg := dom.Sol.Y[eqg]
		plC, _ := dom.Sim.LiqMdl.Calc(z)
		pgC, _ := dom.Sim.GasMdl.Calc(z)
		chk.Scalar(tst, io.Sf("nod %3d : pl(@ %4g)= %6g", nod.Vert.Id, z, pl), 1e-17, pl, plC)
		chk.Scalar(tst, io.Sf("nod %3d : pg(@ %4g)= %6g", nod.Vert.Id, z, pg), 1e-17, pg, pgC)
	}

	// intial values @ integration points
	io.Pforan("initial values @ integration points\n")
	for _, ele := range dom.Elems {
		e := ele.(*seepage.LiquidGas)
		slmax := e.Mdl.Lrm.SlMax()
		for idx, ip := range e.IpsElem {
			s := e.States[idx]
			z := e.Cell.Shp.IpRealCoords(e.X, ip)[1]
			_, ρLC := dom.Sim.LiqMdl.Calc(z)
			_, ρGC := dom.Sim.GasMdl.Calc(z)
			chk.Scalar(tst, io.Sf("sl(@ %18g)= %18g", z, s.A_sl), 1e-17, s.A_sl, slmax)
			chk.Scalar(tst, io.Sf("ρL(@ %18g)= %18g", z, s.A_ρL), 1e-9, s.A_ρL, ρLC)
			chk.Scalar(tst, io.Sf("ρG(@ %18g)= %18g", z, s.A_ρG), 1e-14, s.A_ρG, ρGC)
		}
	}
}

func Test_pp01b(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("pp01b. Liquid-Gas. Run and check Kb")

	// run simulation
	main := fem.NewMain("data/pp01.sim", "", true, false, false, false, chk.Verbose, 0)

	// for debugging Kb
	if true {
		tests.LiquidGas(main, &tests.Kb{
			Tst: tst, Eid: 3, Tol: 1e-6, Verb: chk.Verbose,
			Ni: 1, Nj: 1, ItMin: 1, ItMax: -1, Tmin: 0, Tmax: 2000,
		})
	}

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// plot
	if chk.Verbose {
		//dom := main.Domains[0]
		//ele := dom.Elems[0].(*seepage.LiquidGas)
		//mdl := ele.Mdl
		//Cnd := ele.Mdl.Cnd
		if true {
			//porous.PlotLrm(mdl, "/tmp/gofem", "fig_pp01b_lrm.eps", 30, 101, true, true, true, "", "", "", "", "")
			//conduct.Plot(Cnd, "/tmp/gofem", "fig_pp01b_liq.eps", 101, false, true, true)
		}
		//conduct.Plot(Cnd, "/tmp/gofem", "fig_pp01b_gas.eps", 101, true, true, true)
	}

	// TODO: add check here
}
