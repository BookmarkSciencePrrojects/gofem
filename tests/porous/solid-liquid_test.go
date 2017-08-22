// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"math"
	"sort"
	"testing"

	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/ele/porous"
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/tests"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

func Test_up01a(tst *testing.T) {

	/* this tests simulates seepage flow along a column
	 * by reducing the initial hydrostatic pressure at
	 * at the bottom of the column
	 *
	 *   using mesh from col10m4e2lay.msh
	 *
	 *      Nodes / Tags                       Equations
	 *                              ux uy pl               ux uy pl
	 *     8 o----o----o 9 (-5)     53 54 55  o----o----o  50 51 52
	 *       |   14    |             .  .  .  |  58 59  |   .  .  .
	 *       |  (-1)   |             .  .  .  |         |   .  .  .
	 *    21 o    o    o 22 (-6)    60 61  .  o    o    o  56 57  .
	 *       |   26    |             .  .  .  |  62 63  |   .  .  .
	 *       |         |             .  .  .  |         |   .  .  .
	 *     6 o----o----o 7 (-4)     39 40 41  o----o----o  36 37 38
	 *       |   13    |             .  .  .  |  44 45  |   .  .  .
	 *       |  (-1)   |             .  .  .  |         |   .  .  .
	 *    19 o    o    o 20 (-6)    46 47  .  o    o    o  42 43  .
	 *       |   25    |             .  .  .  |  48 49  |   .  .  .
	 *       |         |             .  .  .  |         |   .  .  .
	 *     4 o----o----o 5 (-3)     25 26 27  o----o----o  22 23 24
	 *       |   12    |             .  .  .  |  30 31  |   .  .  .
	 *       |  (-2)   |             .  .  .  |         |   .  .  .
	 *    17 o    o    o 18 (-6)    32 33  .  o    o    o  28 29  .
	 *       |   24    |             .  .  .  |  34 35  |   .  .  .
	 *       |         |             .  .  .  |         |   .  .  .
	 *     2 o----o----o 3 (-2)      9 10 11  o----o----o   6  7  8
	 *       |   11    |             .  .  .  |  16 17  |   .  .  .
	 *       |  (-2)   |             .  .  .  |         |   .  .  .
	 *    15 o    o    o 16 (-6)    18 19     o    o    o  14 15
	 *       |   23    |             .  .  .  |  20 21  |   .  .  .
	 *       |         |             .  .  .  |         |   .  .  .
	 *     0 o----o----o 1 (-1)      0  1  2  o----o----o   3  4  5
	 *           10                              12 13
	 */

	//tests.Verbose()
	chk.PrintTitle("up01a. Solid-Liquid coupling. Check DOFs")

	// start simulation
	main := fem.NewMain("data/up01.sim", "", true, false, false, false, chk.Verbose, 0)

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

	if true {

		// nodes with pl
		nods_with_pl := map[int]bool{0: true, 2: true, 4: true, 6: true, 8: true, 1: true, 3: true, 5: true, 7: true, 9: true}

		// check dofs
		for _, nod := range dom.Nodes {
			if nods_with_pl[nod.Vert.Id] {
				chk.IntAssert(len(nod.Dofs), 3)
				chk.StrAssert(nod.Dofs[0].Key, "ux")
				chk.StrAssert(nod.Dofs[1].Key, "uy")
				chk.StrAssert(nod.Dofs[2].Key, "pl")
			} else {
				chk.IntAssert(len(nod.Dofs), 2)
				chk.StrAssert(nod.Dofs[0].Key, "ux")
				chk.StrAssert(nod.Dofs[1].Key, "uy")
			}
		}

		// check equations
		nids, eqs := tests.GetNidsEqs(dom)
		chk.Ints(tst, "eqs", eqs, utl.IntRange(10*3+17*2))
		chk.Ints(tst, "nids", nids, []int{
			0, 1, 3, 2, 10, 16, 11, 15, 23,
			5, 4, 18, 12, 17, 24,
			7, 6, 20, 13, 19, 25,
			9, 8, 22, 14, 21, 26,
		})

		// check pmap
		Pmaps := [][]int{
			{2, 5, 8, 11},
			{11, 8, 24, 27},
			{27, 24, 38, 41},
			{41, 38, 52, 55},
		}
		Umaps := [][]int{
			{0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21},
			{9, 10, 6, 7, 22, 23, 25, 26, 16, 17, 28, 29, 30, 31, 32, 33, 34, 35},
			{25, 26, 22, 23, 36, 37, 39, 40, 30, 31, 42, 43, 44, 45, 46, 47, 48, 49},
			{39, 40, 36, 37, 50, 51, 53, 54, 44, 45, 56, 57, 58, 59, 60, 61, 62, 63},
		}
		for i, ele := range dom.Elems {
			e := ele.(*porous.SolidLiquid)
			io.Pfpink("%2d : Pmap = %v\n", e.Id(), e.P.Pmap)
			io.Pfpink("%2d : Umap = %v\n", e.Id(), e.U.Umap)
			chk.Ints(tst, "Pmap", e.P.Pmap, Pmaps[i])
			chk.Ints(tst, "Umap", e.U.Umap, Umaps[i])
		}

		// constraints
		chk.IntAssert(len(dom.EssenBcs.Bcs), 9*2+2+3)
		var ct_ux_eqs []int // equations with ux prescribed [sorted]
		var ct_uy_eqs []int // equations with uy prescribed [sorted]
		var ct_pl_eqs []int // equations with pl prescribed [sorted]
		for _, c := range dom.EssenBcs.Bcs {
			chk.IntAssert(len(c.Eqs), 1)
			eq := c.Eqs[0]
			io.Pfgrey("key=%v eq=%v\n", c.Key, eq)
			switch c.Key {
			case "ux":
				ct_ux_eqs = append(ct_ux_eqs, eq)
			case "uy":
				ct_uy_eqs = append(ct_uy_eqs, eq)
			case "pl":
				ct_pl_eqs = append(ct_pl_eqs, eq)
			default:
				tst.Errorf("key %s is incorrect", c.Key)
			}
		}
		sort.Ints(ct_ux_eqs)
		sort.Ints(ct_uy_eqs)
		sort.Ints(ct_pl_eqs)
		chk.Ints(tst, "equations with ux prescribed", ct_ux_eqs, []int{0, 3, 6, 9, 14, 18, 22, 25, 28, 32, 36, 39, 42, 46, 50, 53, 56, 60})
		chk.Ints(tst, "equations with uy prescribed", ct_uy_eqs, []int{1, 4, 13})
		chk.Ints(tst, "equations with pl prescribed", ct_pl_eqs, []int{2, 5})

	}

	// initial values @ nodes
	io.Pforan("initial values @ nodes\n")
	for _, nod := range dom.Nodes {
		z := nod.Vert.C[1]
		for _, dof := range nod.Dofs {
			u := dom.Sol.Y[dof.Eq]
			switch dof.Key {
			case "ux":
				chk.Float64(tst, io.Sf("nod %3d : ux(@ %4g)= %6g", nod.Vert.Id, z, u), 1e-17, u, 0)
			case "uy":
				chk.Float64(tst, io.Sf("nod %3d : uy(@ %4g)= %6g", nod.Vert.Id, z, u), 1e-17, u, 0)
			case "pl":
				plC, _ := dom.Sim.LiqMdl.Calc(z)
				chk.Float64(tst, io.Sf("nod %3d : pl(@ %4g)= %6g", nod.Vert.Id, z, u), 1e-13, u, plC)
			}
		}
	}

	// layers/column data
	H := dom.Sim.MaxElev
	ν := 0.2
	colTop := fem.ColLayer{
		Zmax: H,
		Zmin: H / 2.0,
		K0:   ν / (1.0 - ν),
		Grav: dom.Sim.Grav0,
		Liq:  dom.Sim.LiqMdl,
		Gas:  dom.Sim.GasMdl,
	}
	colBot := fem.ColLayer{
		Zmax: H / 2.0,
		Zmin: 0.0,
		K0:   ν / (1.0 - ν),
		Grav: dom.Sim.Grav0,
		Liq:  dom.Sim.LiqMdl,
		Gas:  dom.Sim.GasMdl,
	}

	// tolerance for density
	tolRho := 1e-5
	if dom.Sim.Data.NoLBB {
		tolRho = 1e-10
	}

	// intial values @ integration points
	io.Pforan("initial values @ integration points\n")
	for _, ele := range dom.Elems {
		e := ele.(*porous.SolidLiquid)
		slmax := e.P.Mdl.Lrm.SlMax()
		if e.Cell.Tag == -1 {
			colTop.SlMax = slmax
			colTop.Nf0 = e.P.Mdl.Nf0
			colTop.RhoS0 = e.P.Mdl.RhoS0
		} else {
			colBot.SlMax = slmax
			colBot.Nf0 = e.P.Mdl.Nf0
			colBot.RhoS0 = e.P.Mdl.RhoS0
		}
		for idx, ip := range e.U.IpsElem {
			s := e.P.States[idx]
			z := e.P.Cell.Shp.IpRealCoords(e.P.X, ip)[1]
			_, ρLC := dom.Sim.LiqMdl.Calc(z)
			chk.AnaNum(tst, io.Sf("sl(@%6.4f)", z), 1e-17, s.A_sl, slmax, chk.Verbose)
			chk.AnaNum(tst, io.Sf("ρL(@%6.4f)", z), tolRho, s.A_ρL, ρLC, chk.Verbose)
		}
	}

	// vertical stress @ zmax of bottom layer
	_, _, _, _, _, colBot.SigV = colTop.Calc(colTop.Zmin)

	// check stresses
	io.Pforan("initial stresses @ integration points\n")
	for _, ele := range dom.Elems {
		e := ele.(*porous.SolidLiquid)
		col := &colTop
		if e.Cell.Tag == -2 {
			col = &colBot
		}
		for idx, ip := range e.U.IpsElem {
			z := e.U.Cell.Shp.IpRealCoords(e.U.X, ip)[1]
			σe := e.U.States[idx].Sig
			pl, _, _, _, _, σV := col.Calc(z)
			sl := col.SlMax
			p := sl * pl
			sve := -σV + p
			she := sve * col.K0
			if math.Abs(σe[2]-σe[0]) > 1e-17 {
				tst.Errorf("[1;31mσx is not equal to σz: %g != %g[0m\n", σe[2], σe[0])
				return
			}
			if math.Abs(σe[3]) > 1e-17 {
				tst.Errorf("[1;31mσxy is not equal to zero: %g != 0[0m\n", σe[3])
				return
			}
			chk.AnaNum(tst, io.Sf("sy(z=%11.8f)", z), 1e-4, σe[1], sve, chk.Verbose)
			chk.AnaNum(tst, io.Sf("sx(z=%11.8f)", z), 1e-4, σe[0], she, chk.Verbose)
		}
	}
}

func Test_up01b(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("up01b. Solid-Liquid coupling. Run")

	// start simulation
	main := fem.NewMain("data/up01.sim", "", true, true, false, false, chk.Verbose, 0)

	// for debugging Kb
	if true {
		tests.SolidLiquid(main, &tests.Kb{
			Tst: tst, Eid: 3, Tol: 1e-8, Tol2: 1e-5, Verb: chk.Verbose,
			Ni: 1, Nj: 1, ItMin: 1, ItMax: -1, Tmin: -1, Tmax: -1,
		})
	}

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
	}

	// (low-level) post-processing
	if chk.Verbose {

		// start prost-processing structure
		postproc := fem.NewMain("data/up01.sim", "", false, false, true, false, chk.Verbose, 0)

		// set stage
		err := postproc.SetStage(0)
		if err != nil {
			chk.Panic("cannot set stage:\n%v", err)
		}

		// initialise solution vectors
		err = postproc.ZeroStage(0, true)
		if err != nil {
			chk.Panic("cannot initialise solution vectors:\n%v", err)
		}

		// domain and summary
		dom := postproc.Domains[0]
		sum := postproc.Summary
		e := dom.Elems[3].(*porous.SolidLiquid)

		// for each output time
		nip := len(e.OutIpCoords())
		ntout := len(sum.OutTimes)
		Pc, Sl := utl.Alloc(nip, ntout), utl.Alloc(nip, ntout)
		for tidx, _ := range sum.OutTimes {

			// input results into domain
			err := dom.Read(sum, tidx, 0, true)
			if err != nil {
				chk.Panic("cannot load results into domain:\n%v", err)
			}

			// results @ ip
			M := ele.NewIpsMap()
			e.OutIpVals(M, dom.Sol)
			for i := 0; i < nip; i++ {
				Pc[i][tidx], Sl[i][tidx] = -M.Get("pl", i), M.Get("sl", i)
			}
		}

		// plot
		//Lrm := e.P.Mdl.Lrm
		//retention.Plot(Lrm, 0, Lrm.SlMax(), 30, 101, false, "'k-^', markerfacecolor='white', ms=5, markevery=10", "", "model")
		//M := []string{".", "+", "x", "1", "2", "3", "4", "|", "_"}
		//for i := 0; i < nip; i++ {
		//plt.Plot(Pc[i], Sl[i], io.Sf("'b.', clip_on=0, markevery=%d, marker='%s'", 3+i, M[i]))
		//}
		//plt.Gll("$p_c$", "$s_{\\ell}$", "")
		//plt.SaveD("/tmp/gofem", "fig_up01.eps")
	}
}
