// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"sort"
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func Test_bh16a(tst *testing.T) {

	/*  solid bracket with thickness = 0.25
	 *
	 *          1     -10                connectivity:
	 *   (-100) o'-,__                    eid :  verts
	 *          |     '-,__ 3   -10         0 : 0, 2, 3
	 *          |        ,'o-,__            1 : 3, 1, 0
	 *          |  1   ,'  |    '-,__ 5     2 : 2, 4, 5
	 *          |    ,'    |  3   ,-'o      3 : 5, 3, 2
	 *          |  ,'  0   |   ,-'   |
	 *          |,'        |,-'   2  |   constraints:
	 *   (-100) o----------o---------o    -100 : fixed on x and y
	 *          0          2         4
	 */

	//verbose()
	chk.PrintTitle("bh16a. bracket. check DOFs")

	// start simulation
	analysis := NewFEM("data/bh16.sim", "", true, false, false, false, chk.Verbose, 0)

	// set stage
	err := analysis.SetStage(0)
	if err != nil {
		tst.Errorf("SetStage failed:\n%v", err)
		return
	}

	// initialise solution vectros
	err = analysis.ZeroStage(0, true)
	if err != nil {
		tst.Errorf("ZeroStage failed:\n%v", err)
		return
	}

	// domain
	dom := analysis.Domains[0]
	io.Pforan("dom.elems = %v\n", dom.Elems)

	// nodes and elements
	chk.IntAssert(len(dom.Nodes), 6)
	chk.IntAssert(len(dom.Elems), 4)

	// check dofs
	for _, nod := range dom.Nodes {
		chk.IntAssert(len(nod.Dofs), 2)
	}

	// check equations
	nids, eqs := get_nids_eqs(dom)
	chk.Ints(tst, "nids", nids, []int{0, 2, 3, 1, 4, 5})
	chk.Ints(tst, "eqs", eqs, []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11})

	// check solution arrays
	ny := 6 * 2
	nλ := 4
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
		{2, 3, 8, 9, 10, 11},
		{10, 11, 4, 5, 2, 3},
	}
	for i, ele := range dom.Elems {
		e := ele.(*ElemU)
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
	chk.Ints(tst, "constrained ux equations", ct_ux_eqs, []int{0, 6})
	chk.Ints(tst, "constrained uy equations", ct_uy_eqs, []int{1, 7})

	// check ip data
	for _, ele := range dom.Elems {
		e := ele.(*ElemU)
		d := e.OutIpsData()
		chk.IntAssert(len(d), 1)
		vals := d[0].Calc(dom.Sol)
		chk.IntAssert(len(vals), 4)
		for key, val := range vals {
			io.Pfyel("key=%v => val=%v\n", key, val)
		}
	}
}

func Test_bh16b(tst *testing.T) {

	//verbose()
	chk.PrintTitle("bh16b. bracket. run")

	// start simulation
	analysis := NewFEM("data/bh16.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 1e-12
	tolu := 1e-15
	tols := 1e-12
	TestingCompareResultsU(tst, "data/bh16.sim", "cmp/bh16.cmp", "", tolK, tolu, tols, skipK, chk.Verbose)
}

func Test_bh14a(tst *testing.T) {

	//verbose()
	chk.PrintTitle("bh14a. using RunAll")

	// start simulation
	analysis := NewFEM("data/bh14.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 1e-10
	tolu := 1e-15
	tols := 1e-13
	TestingCompareResultsU(tst, "data/bh14.sim", "cmp/bh14.cmp", "", tolK, tolu, tols, skipK, chk.Verbose)
}

func Test_bh14b(tst *testing.T) {

	//verbose()
	chk.PrintTitle("bh14b. truss. using SolveOneStage")

	// start simulation
	analysis := NewFEM("data/bh14.sim", "", true, true, false, false, chk.Verbose, 0)

	// set stage
	err := analysis.SetStage(0)
	if err != nil {
		tst.Errorf("SetStage failed:\n%v", err)
		return
	}

	// run
	err = analysis.SolveOneStage(0, true)
	if err != nil {
		tst.Error("SolveOneStage failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 1e-10
	tolu := 1e-15
	tols := 1e-13
	TestingCompareResultsU(tst, "data/bh14.sim", "cmp/bh14.cmp", "", tolK, tolu, tols, skipK, chk.Verbose)
}

func Test_bh14c(tst *testing.T) {

	//verbose()
	chk.PrintTitle("bh14c. truss. using go-routines")

	// channels
	nch := 3
	done := make(chan int, nch)

	// allocate structures and set stage
	analyses := make([]*FEM, nch)
	for i := 0; i < nch; i++ {
		analyses[i] = NewFEM("data/bh14.sim", io.Sf("ch%d", i), true, true, false, false, false, i)
		err := analyses[i].SetStage(0)
		if err != nil {
			tst.Errorf("SetStage failed:\n%v", err)
			return
		}
	}

	// run all analyses
	for i := 0; i < nch; i++ {
		go func(analysis *FEM) {
			err := analysis.SolveOneStage(0, true)
			if err != nil {
				io.Pfred("SolveOneStage failed:\n%v", err)
			}
			done <- 1
		}(analyses[i])
	}

	// wait
	for i := 0; i < nch; i++ {
		<-done
	}

	// check
	skipK := false
	tolK := 1e-10
	tolu := 1e-15
	tols := 1e-13
	for i := 0; i < nch; i++ {
		TestingCompareResultsU(tst, "data/bh14.sim", "cmp/bh14.cmp", io.Sf("ch%d", i), tolK, tolu, tols, skipK, chk.Verbose)
	}
}

func Test_bh14erod(tst *testing.T) {

	//verbose()
	chk.PrintTitle("bh14erod. truss. using ElasticRod")

	// start simulation
	analysis := NewFEM("data/bh14erod.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 1e-10
	tolu := 1e-15
	tols := 1e-13
	TestingCompareResultsU(tst, "data/bh14erod.sim", "cmp/bh14.cmp", "", tolK, tolu, tols, skipK, chk.Verbose)
}
