// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"math"
	"sort"
	"testing"

	"github.com/cpmech/gofem/ele/porous"
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/tests"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

func Test_upp01a(tst *testing.T) {

	/* this tests simulates seepage flow along a column
	 * by reducing the initial hydrostatic pressure at
	 * at the bottom of the column
	 *
	 *   using mesh from col10m4e2lay.msh
	 *
	 *      Nodes / Tags                         Equations
	 *                              ux uy pl pg              ux uy pl pg
	 *     8 o----o----o 9 (-5)     62 63 64 65 o----o----o  58 59 60 61
	 *       |   14    |             .  .  .  . |  68 69  |   .  .  .  .
	 *       |  (-1)   |             .  .  .  . |         |   .  .  .  .
	 *    21 o    o    o 22 (-6)    70 71  .  . o    o    o  66 67  .  .
	 *       |   26    |             .  .  .  . |  72 73  |   .  .  .  .
	 *       |         |             .  .  .  . |         |   .  .  .  .
	 *     6 o----o----o 7 (-4)     46 47 48 49 o----o----o  42 43 44 45
	 *       |   13    |             .  .  .  . |  52 53  |   .  .  .  .
	 *       |  (-1)   |             .  .  .  . |         |   .  .  .  .
	 *    19 o    o    o 20 (-6)    54 55  .  . o    o    o  50 51  .  .
	 *       |   25    |             .  .  .  . |  56 57  |   .  .  .  .
	 *       |         |             .  .  .  . |         |   .  .  .  .
	 *     4 o----o----o 5 (-3)     30 31 32 33 o----o----o  26 27 28 29
	 *       |   12    |             .  .  .  . |  36 37  |   .  .  .  .
	 *       |  (-2)   |             .  .  .  . |         |   .  .  .  .
	 *    17 o    o    o 18 (-6)    38 39  .  . o    o    o  34 35  .  .
	 *       |   24    |             .  .  .  . |  40 41  |   .  .  .  .
	 *       |         |             .  .  .  . |         |   .  .  .  .
	 *     2 o----o----o 3 (-2)     12 13 14 15 o----o----o   8  9 10 11
	 *       |   11    |             .  .  .  . |  20 21  |   .  .  .  .
	 *       |  (-2)   |             .  .  .  . |         |   .  .  .  .
	 *    15 o    o    o 16 (-6)    22 23     . o    o    o  18 19     .
	 *       |   23    |             .  .  .  . |  24 25  |   .  .  .  .
	 *       |         |             .  .  .  . |         |   .  .  .  .
	 *     0 o----o----o 1 (-1)      0  1  2  3 o----o----o   4  5  6  7
	 *           10                                16 17
	 */

	//tests.Verbose()
	chk.PrintTitle("upp01a. Solid-Liquid-Gas coupling. Check DOFs and BCs")

	// start simulation
	main := fem.NewMain("data/upp01.sim", "", true, false, false, false, chk.Verbose, 0)

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

		// nodes with pl or pg
		nods_with_p := map[int]bool{0: true, 2: true, 4: true, 6: true, 8: true, 1: true, 3: true, 5: true, 7: true, 9: true}

		// check dofs
		for _, nod := range dom.Nodes {
			if nods_with_p[nod.Vert.Id] {
				chk.IntAssert(len(nod.Dofs), 4)
				chk.StrAssert(nod.Dofs[0].Key, "ux")
				chk.StrAssert(nod.Dofs[1].Key, "uy")
				chk.StrAssert(nod.Dofs[2].Key, "pl")
				chk.StrAssert(nod.Dofs[3].Key, "pg")
			} else {
				chk.IntAssert(len(nod.Dofs), 2)
				chk.StrAssert(nod.Dofs[0].Key, "ux")
				chk.StrAssert(nod.Dofs[1].Key, "uy")
			}
		}

		// check equations
		nids, eqs := tests.GetNidsEqs(dom)
		chk.Ints(tst, "eqs", eqs, utl.IntRange(10*4+17*2))
		chk.Ints(tst, "nids", nids, []int{
			0, 1, 3, 2, 10, 16, 11, 15, 23,
			5, 4, 18, 12, 17, 24,
			7, 6, 20, 13, 19, 25,
			9, 8, 22, 14, 21, 26,
		})

		// check maps
		Plmaps := [][]int{
			{2, 6, 10, 14},
			{14, 10, 28, 32},
			{32, 28, 44, 48},
			{48, 44, 60, 64},
		}
		Pgmaps := [][]int{
			{3, 7, 11, 15},
			{15, 11, 29, 33},
			{33, 29, 45, 49},
			{49, 45, 61, 65},
		}
		Umaps := [][]int{
			{0, 1, 4, 5, 8, 9, 12, 13, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25},
			{12, 13, 8, 9, 26, 27, 30, 31, 20, 21, 34, 35, 36, 37, 38, 39, 40, 41},
			{30, 31, 26, 27, 42, 43, 46, 47, 36, 37, 50, 51, 52, 53, 54, 55, 56, 57},
			{46, 47, 42, 43, 58, 59, 62, 63, 52, 53, 66, 67, 68, 69, 70, 71, 72, 73},
		}
		for i, ele := range dom.Elems {
			e := ele.(*porous.SolidLiquidGas)
			chk.Ints(tst, "Plmap", e.P.Plmap, Plmaps[i])
			chk.Ints(tst, "Pgmap", e.P.Pgmap, Pgmaps[i])
			chk.Ints(tst, "Umap ", e.U.Umap, Umaps[i])
		}

		// constraints
		chk.IntAssert(len(dom.EssenBcs.Bcs), 9*2+3+4+2)
		var ct_ux_eqs []int // equations with ux prescribed [sorted]
		var ct_uy_eqs []int // equations with uy prescribed [sorted]
		var ct_pl_eqs []int // equations with pl prescribed [sorted]
		var ct_pg_eqs []int // equations with pg prescribed [sorted]
		for _, c := range dom.EssenBcs.Bcs {
			chk.IntAssert(len(c.Eqs), 1)
			eq := c.Eqs[0]
			switch c.Key {
			case "ux":
				ct_ux_eqs = append(ct_ux_eqs, eq)
			case "uy":
				ct_uy_eqs = append(ct_uy_eqs, eq)
			case "pl":
				ct_pl_eqs = append(ct_pl_eqs, eq)
			case "pg":
				ct_pg_eqs = append(ct_pg_eqs, eq)
			default:
				tst.Errorf("key %s is incorrect", c.Key)
			}
		}
		sort.Ints(ct_ux_eqs)
		sort.Ints(ct_uy_eqs)
		sort.Ints(ct_pl_eqs)
		sort.Ints(ct_pg_eqs)
		chk.Ints(tst, "equations with ux prescribed", ct_ux_eqs, []int{0, 4, 8, 12, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62, 66, 70})
		chk.Ints(tst, "equations with uy prescribed", ct_uy_eqs, []int{1, 5, 17})
		chk.Ints(tst, "equations with pl prescribed", ct_pl_eqs, []int{2, 6})
		chk.Ints(tst, "equations with pg prescribed", ct_pg_eqs, []int{3, 7, 61, 65})

	}

	// initial values @ nodes
	io.Pforan("initial values @ nodes\n")
	for _, nod := range dom.Nodes {
		z := nod.Vert.C[1]
		for _, dof := range nod.Dofs {
			u := dom.Sol.Y[dof.Eq]
			switch dof.Key {
			case "ux":
				chk.Scalar(tst, io.Sf("nod %3d : ux(@ %4g)= %6g", nod.Vert.Id, z, u), 1e-17, u, 0)
			case "uy":
				chk.Scalar(tst, io.Sf("nod %3d : uy(@ %4g)= %6g", nod.Vert.Id, z, u), 1e-17, u, 0)
			case "pl":
				plC, _ := dom.Sim.LiqMdl.Calc(z)
				chk.Scalar(tst, io.Sf("nod %3d : pl(@ %4g)= %6g", nod.Vert.Id, z, u), 1e-13, u, plC)
			case "pg":
				pgC, _ := dom.Sim.GasMdl.Calc(z)
				chk.Scalar(tst, io.Sf("nod %3d : pg(@ %4g)= %6g", nod.Vert.Id, z, u), 1e-13, u, pgC)
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
		e := ele.(*porous.SolidLiquidGas)
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
			_, ρGC := dom.Sim.GasMdl.Calc(z)
			chk.AnaNum(tst, io.Sf("sl(@%6.4f)", z), 1e-17, s.A_sl, slmax, chk.Verbose)
			chk.AnaNum(tst, io.Sf("ρL(@%6.4f)", z), tolRho, s.A_ρL, ρLC, chk.Verbose)
			chk.AnaNum(tst, io.Sf("ρG(@%6.4f)", z), tolRho, s.A_ρG, ρGC, chk.Verbose)
		}
	}

	// vertical stress @ zmax of bottom layer
	_, _, _, _, _, colBot.SigV = colTop.Calc(colTop.Zmin)

	// check stresses
	io.Pforan("initial stresses @ integration points\n")
	for _, ele := range dom.Elems {
		e := ele.(*porous.SolidLiquidGas)
		col := &colTop
		if e.Cell.Tag == -2 {
			col = &colBot
		}
		for idx, ip := range e.U.IpsElem {
			z := e.U.Cell.Shp.IpRealCoords(e.U.X, ip)[1]
			σe := e.U.States[idx].Sig
			pl, pg, _, _, _, σV := col.Calc(z)
			sl := col.SlMax
			sg := 1.0 - sl
			p := sl*pl + sg*pg
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

func Test_upp01b(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("upp01b. Solid-Liquid-Gas coupling. Run and check Kb")

	// start simulation
	main := fem.NewMain("data/upp01.sim", "", true, false, false, false, chk.Verbose, 0)

	// for debugging Kb
	if true {
		tests.SolidLiquidGas(main, &tests.Kb{
			Tst: tst, Eid: 3, Tol: 1e-7, Tol2: 1e-5, Verb: chk.Verbose,
			Ni: 1, Nj: 1, ItMin: 1, ItMax: -1, Tmin: -1, Tmax: -1,
		})
	}

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
	}
}
