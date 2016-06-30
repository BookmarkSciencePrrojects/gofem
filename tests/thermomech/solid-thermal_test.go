// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"sort"
	"testing"
	"github.com/cpmech/gofem/ele/thermomech"
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/tests"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

func Test_ut01a(tst *testing.T) {

	/*Test thermomechanical coupling
	 *
	 *   using mesh from col10m4e2lay.msh
	 *
	 *      Nodes / Tags                       Equations
	 *                              ux uy  t               ux uy  t
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

	tests.Verbose()
	chk.PrintTitle("ut01a")

	// start simulation
	main := fem.NewMain("data/ut01.sim", "", true, false, false, false, chk.Verbose, 0)

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
		nods_with_t := map[int]bool{0: true, 2: true, 4: true, 6: true, 8: true, 1: true, 3: true, 5: true, 7: true, 9: true}

		// check dofs
		for _, nod := range dom.Nodes {
			if nods_with_t[nod.Vert.Id] {
				chk.IntAssert(len(nod.Dofs), 3)
				chk.StrAssert(nod.Dofs[0].Key, "ux")
				chk.StrAssert(nod.Dofs[1].Key, "uy")
				chk.StrAssert(nod.Dofs[2].Key, "temp")
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
		Tmaps := [][]int{
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
			e := ele.(*thermomech.SolidThermal)
			io.Pfpink("%2d : Tmap = %v\n", e.Id(), e.Tmap)
			io.Pfpink("%2d : Umap = %v\n", e.Id(), e.Umap)
			chk.Ints(tst, "Tmap", e.Tmap, Tmaps[i])
			chk.Ints(tst, "Umap", e.Umap, Umaps[i])
		}

		// constraints
		//chk.IntAssert(len(dom.EssenBcs.Bcs), 3+3+5+5)
		var ct_ux_eqs []int // equations with ux prescribed [sorted]
		var ct_uy_eqs []int // equations with uy prescribed [sorted]
		var ct_t_eqs []int // equations with t prescribed [sorted]
		for _, c := range dom.EssenBcs.Bcs {
			chk.IntAssert(len(c.Eqs), 1)
			eq := c.Eqs[0]
			io.Pfgrey("key=%v eq=%v\n", c.Key, eq)
			switch c.Key {
			case "ux":
				ct_ux_eqs = append(ct_ux_eqs, eq)
			case "uy":
				ct_uy_eqs = append(ct_uy_eqs, eq)
			case "temp":
				ct_t_eqs = append(ct_t_eqs, eq)
			default:
				tst.Errorf("key %s is incorrect", c.Key)
			}
		}
		sort.Ints(ct_ux_eqs)
		sort.Ints(ct_uy_eqs)
		sort.Ints(ct_t_eqs)
		chk.Ints(tst, "equations with ux prescribed", ct_ux_eqs, []int{0, 3, 12})
		chk.Ints(tst, "equations with uy prescribed", ct_uy_eqs, []int{1, 4, 13})
		chk.Ints(tst, "equations with t prescribed", ct_t_eqs, []int{2, 5, 8, 11, 24, 27, 38, 41, 52, 55})

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
			case "temp":
				chk.Scalar(tst, io.Sf("nod %3d : t(@ %4g)= %6g", nod.Vert.Id, z, u), 1e-17, u, 0)
			}
		}
	}
}

func Test_ut01b(tst *testing.T) {

	tests.Verbose()
	chk.PrintTitle("ut01b")

	// start simulation
	main := fem.NewMain("data/ut01.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
	}

	// domain and element
	dom := main.Domains[0]

	// check displacements and temperatures
	for _, n := range dom.Nodes {
		eqx := n.GetEq("ux")
		eqy := n.GetEq("uy")
		eqt := n.GetEq("temp")
		//fmt.Printf("id=%v, dofs=%v", n, len(n.Dofs))
		if len(n.Dofs) == 2 {
			u := []float64{dom.Sol.Y[eqx], dom.Sol.Y[eqy]}
			io.Pfyel("x=%4.3v   u=%10.4v\n", n.Vert.C, u)
		} else {
			u := []float64{dom.Sol.Y[eqx], dom.Sol.Y[eqy], dom.Sol.Y[eqt]}
			io.Pfyel("x=%4.3v   u=%10.4v\n", n.Vert.C, u)
		}
	}
	//Analytical solution from Richard B. Hetnarski, Thermal Stresses - Advanced Theory and Applications, p.228
	io.Pfyel("\nANALYTICAL RESULT u[0]:\nx=[0   10]   u= %10.4f\n", 1.0/1000*(100.0-0.0)*10.0*10.0/2.0/2.5)
}

