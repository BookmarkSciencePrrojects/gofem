// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"sort"
	"testing"

	"github.com/cpmech/gofem/ele/solid"
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/tests"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func Test_sgm52a(tst *testing.T) {

	/* Smith, Griffiths and Margetts (5th ed) Figure 5.2 p173
	 *
	 *          0.25       0.5      0.25 kN/m
	 *            ↓         ↓         ↓
	 *    ---    ▷0---------1---------2
	 *     |      |       ,'|       ,'|   E = 1e6 kN/m²
	 *     |      |  0  ,'  |  2  ,'  |   ν = 0.3
	 *     |      |   ,'    |   ,'    |
	 *            | ,'   1  | ,'  3   |   connectivity:
	 *    1 m    ▷3'--------4'--------5     0 : 1 0 3
	 *            |       ,'|       ,'|     1 : 3 4 1
	 *     |      |  4  ,'  |  6  ,'  |     2 : 2 1 4
	 *     |      |   ,'    |   ,'    |     3 : 4 5 2
	 *     |      | ,'   5  | ,'   7  |     4 : 4 3 6
	 *    ---    ▷6'--------7'--------8     5 : 6 7 4
	 *            △         △         △     6 : 5 4 7
	 *                                      7 : 7 8 5
	 *            |------- 1 m -------|
	 */

	//tests.Verbose()
	chk.PrintTitle("sgm52a. plane strain tri3. check DOFs")

	// start simulation
	main := fem.NewMain("data/sgm52.sim", "", true, false, false, false, chk.Verbose, 0)

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
	chk.IntAssert(len(dom.Nodes), 9)
	chk.IntAssert(len(dom.Elems), 8)

	// check dofs
	for _, nod := range dom.Nodes {
		chk.IntAssert(len(nod.Dofs), 2)
	}

	// check equations
	nids, eqs := tests.GetNidsEqs(dom)
	chk.Ints(tst, "nids", nids, []int{1, 0, 3, 4, 2, 5, 6, 7, 8})
	chk.Ints(tst, "eqs", eqs, []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17})

	// check solution arrays
	ny := 9 * 2
	nλ := 6
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
	umaps := [][]int{
		{0, 1, 2, 3, 4, 5},
		{4, 5, 6, 7, 0, 1},
		{8, 9, 0, 1, 6, 7},
		{6, 7, 10, 11, 8, 9},
		{6, 7, 4, 5, 12, 13},
		{12, 13, 14, 15, 6, 7},
		{10, 11, 6, 7, 14, 15},
		{14, 15, 16, 17, 10, 11},
	}
	for i, ele := range dom.Elems {
		e := ele.(*solid.Solid)
		io.Pforan("e%d.umap = %v\n", e.Id(), e.Umap)
		chk.Ints(tst, "umap", e.Umap, umaps[i])
	}

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
		}
	}
	sort.Ints(ct_ux_eqs)
	sort.Ints(ct_uy_eqs)
	chk.Ints(tst, "constrained ux equations", ct_ux_eqs, []int{2, 4, 12})
	chk.Ints(tst, "constrained uy equations", ct_uy_eqs, []int{13, 15, 17})

	// point loads
	chk.IntAssert(len(dom.PtNatBcs.Bcs), 3)
	chk.StrAssert(dom.PtNatBcs.Bcs[0].Key, "fuy")
	chk.StrAssert(dom.PtNatBcs.Bcs[1].Key, "fuy")
	chk.StrAssert(dom.PtNatBcs.Bcs[2].Key, "fuy")
	chk.IntAssert(dom.PtNatBcs.Bcs[0].Eq, 3)
	chk.IntAssert(dom.PtNatBcs.Bcs[1].Eq, 1)
	chk.IntAssert(dom.PtNatBcs.Bcs[2].Eq, 9)
}

func Test_sgm52b(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("sgm52b. plane strain tri3. run")

	// run simulation
	main := fem.NewMain("data/sgm52.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 1e-9
	tolu := 1e-17
	tols := 1.56e-15
	tests.CompareResults(tst, "data/sgm52.sim", "cmp/sgm52.cmp", "", tolK, tolu, tols, skipK, chk.Verbose, nil)
}

func Test_sgm57(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("sgm57. plane strain tri15. qn given")

	// run simulation
	main := fem.NewMain("data/sgm57.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 0.35
	tolu := 2e-9
	tols := 0.0002
	tests.CompareResults(tst, "data/sgm57.sim", "cmp/sgm57.cmp", "", tolK, tolu, tols, skipK, chk.Verbose, nil)
}

func Test_sgm511(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("sgm511. plane strain qua4. disp given")

	// run simulation
	main := fem.NewMain("data/sgm511.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 0.1
	tolu := 3e-14
	tols := 1.56e-7
	tests.CompareResults(tst, "data/sgm511.sim", "cmp/sgm511.cmp", "", tolK, tolu, tols, skipK, chk.Verbose, nil)
}

func Test_sgm515(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("sgm515. plane strain qua8. qn given")

	// run simulation
	main := fem.NewMain("data/sgm515.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// displacement @ top
	dom := main.Domains[0]
	nod := dom.Vid2node[0]
	eqy := nod.GetEq("uy")
	uy := dom.Sol.Y[eqy]
	uy_cor := -5.310671749739340E-06
	io.Pforan("uy @ top = %v (%v)\n", uy, uy_cor)
	chk.Scalar(tst, "uy @ top", 1e-12, uy, uy_cor)

	// check
	if true {
		skipK := false
		tolK := 0.15
		tolu := 3e-13
		tols := 3e-8
		tests.CompareResults(tst, "data/sgm515.sim", "cmp/sgm515.cmp", "", tolK, tolu, tols, skipK, chk.Verbose, nil)
	}
}

func Test_sgm527(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("sgm527. plane strain qua9. qn given")

	// run simulation
	main := fem.NewMain("data/sgm527.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// displacement @ top
	dom := main.Domains[0]
	nod := dom.Vid2node[0]
	eqy := nod.GetEq("uy")
	uy := dom.Sol.Y[eqy]
	uy_cor := -5.298917281673461E-06
	io.Pforan("uy @ top = %v (%v)\n", uy, uy_cor)
	chk.Scalar(tst, "uy @ top", 1e-12, uy, uy_cor)

	// check
	if true {
		skipK := true
		tolK := 1e-17
		tolu := 3e-13
		tols := 5e-8
		tests.CompareResults(tst, "data/sgm527.sim", "cmp/sgm527.cmp", "", tolK, tolu, tols, skipK, chk.Verbose, nil)
	}
}

func Test_sgm517(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("sgm517. axisymmetric qua4")

	// run simulation
	main := fem.NewMain("data/sgm517.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 0.0036
	tolu := 1e-6
	tols := 1e-4
	tests.CompareResults(tst, "data/sgm517.sim", "cmp/sgm517.cmp", "", tolK, tolu, tols, skipK, chk.Verbose, nil)
}

func Test_sgm524(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("sgm524. 3d hex8")

	// run simulation
	main := fem.NewMain("data/sgm524.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := true
	tolK := 1e-17
	tolu := 1e-8
	tols := 1e-7
	tests.CompareResults(tst, "data/sgm524.sim", "cmp/sgm524.cmp", "", tolK, tolu, tols, skipK, chk.Verbose, nil)
}

func Test_sgm530(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("sgm530. 3d tet4")

	// run simulation
	main := fem.NewMain("data/sgm530.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := true
	tolK := 1e-17
	tolu := 1e-17
	tols := 1e-15
	tests.CompareResults(tst, "data/sgm530.sim", "cmp/sgm530.cmp", "", tolK, tolu, tols, skipK, chk.Verbose, nil)
}

func Test_sgm422(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("sgm422. small 3D frame")

	// run simulation
	main := fem.NewMain("data/sgm422.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	dom := main.Domains[0]
	nod := dom.Nodes[1]
	eqy := nod.GetEq("uy")
	uy := dom.Sol.Y[eqy]
	io.Pforan("uy @ nod 2 = %v\n", uy)

	// check
	skipK := true
	tolK := 1e-17
	tolu := 1e-17
	tols := 1e-15
	tests.CompareResults(tst, "data/sgm422.sim", "cmp/sgm422.cmp", "", tolK, tolu, tols, skipK, chk.Verbose, nil)
}
