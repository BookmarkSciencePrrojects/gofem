// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"math"
	"sort"
	"testing"

	"github.com/cpmech/gofem/ele/solid"
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/tests"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

func Test_spo751a(tst *testing.T) {

	/* de Souza Neto, Perić and Owen, ex 7.5.1 p244
	 *
	 *                       22
	 *                        .
	 *                  19  ,' `.
	 *                    ,'     '.
	 *              17  ,'         \
	 *                .'            \
	 *           14 ,' `.            \ 21
	 *         12 ,'     \            '
	 *       9  .'        \            '
	 *     7  ,' `.        \ 16         '
	 *   4  .'     \        .           `
	 *  2  ' `.     \ 11     .          |
	 *    `.   \ 6   .       |          |
	 *     1.   .    |       |          |
	 *      |   |    |       |          |
	 *      -----------------------------
	 *      0 3 5 8 10  13  15    18   20
	 */

	//tests.Verbose()
	chk.PrintTitle("spo751a. cylinder expansion. check DOFs")

	// start simulation
	main := fem.NewMain("data/spo751.sim", "", true, false, false, false, chk.Verbose, 0)

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
	chk.IntAssert(len(dom.Nodes), 23)
	chk.IntAssert(len(dom.Elems), 4)

	// check dofs
	for _, nod := range dom.Nodes {
		chk.IntAssert(len(nod.Dofs), 2)
	}

	// check equations
	nids, eqs := tests.GetNidsEqs(dom)
	chk.Ints(tst, "nids", nids, []int{0, 5, 7, 2, 3, 6, 4, 1, 10, 12, 8, 11, 9, 15, 17, 13, 16, 14, 20, 22, 18, 21, 19})
	chk.Ints(tst, "eqs", eqs, utl.IntRange(23*2))

	// check solution arrays
	ny := 23 * 2
	nλ := 9 + 9
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
		{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
		{2, 3, 16, 17, 18, 19, 4, 5, 20, 21, 22, 23, 24, 25, 10, 11},
		{16, 17, 26, 27, 28, 29, 18, 19, 30, 31, 32, 33, 34, 35, 22, 23},
		{26, 27, 36, 37, 38, 39, 28, 29, 40, 41, 42, 43, 44, 45, 32, 33},
	}
	for i, ele := range dom.Elems {
		e := ele.(*solid.Solid)
		io.Pforan("e%d.umap = %v\n", e.Id(), e.Umap)
		chk.Ints(tst, "umap", e.Umap, umaps[i])
	}

	// constraints
	chk.IntAssert(len(dom.EssenBcs.Bcs), nλ)
	var ct_uy_eqs []int // constrained uy equations [sorted]
	var ct_incsup_xeqs []int
	var ct_incsup_yeqs []int
	αrad := 120.0 * math.Pi / 180.0
	cα, sα := math.Cos(αrad), math.Sin(αrad)
	for _, c := range dom.EssenBcs.Bcs {
		chk.IntAssertLessThanOrEqualTo(1, len(c.Eqs)) // 1 ≤ neqs
		io.Pforan("c.Key=%s c.Eqs=%v\n", c.Key, c.Eqs)
		if len(c.Eqs) == 1 {
			if c.Key == "uy" {
				ct_uy_eqs = append(ct_uy_eqs, c.Eqs[0])
				continue
			}
		} else {
			if c.Key == "incsup" {
				ct_incsup_xeqs = append(ct_incsup_xeqs, c.Eqs[0])
				ct_incsup_yeqs = append(ct_incsup_yeqs, c.Eqs[1])
				chk.AnaNum(tst, "cos(α)", 1e-15, c.ValsA[0], cα, false)
				chk.AnaNum(tst, "sin(α)", 1e-15, c.ValsA[1], sα, false)
				continue
			}
		}
		tst.Errorf("key %s is incorrect", c.Key)
	}
	sort.Ints(ct_uy_eqs)
	sort.Ints(ct_incsup_xeqs)
	sort.Ints(ct_incsup_yeqs)
	chk.Ints(tst, "constrained uy equations", ct_uy_eqs, []int{1, 3, 9, 17, 21, 27, 31, 37, 41})
	chk.Ints(tst, "incsup x equations", ct_incsup_xeqs, []int{4, 6, 12, 18, 24, 28, 34, 38, 44})
	chk.Ints(tst, "incsup y equations", ct_incsup_yeqs, []int{5, 7, 13, 19, 25, 29, 35, 39, 45})
}

func Test_spo751b(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("spo751b. cylinder expansion. run")

	// run simulation
	main := fem.NewMain("data/spo751.sim", "", true, true, false, false, chk.Verbose, 0)

	// for debugging Kb
	if true {
		tests.Solid(main, &tests.Kb{
			Tst: tst, Eid: 3, Tol: 1e-5, Verb: chk.Verbose,
			Ni: 1, Nj: 1, ItMin: 1, ItMax: -1, Tmin: 0.89, Tmax: 0.96,
		})
	}

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	if true {
		verb := false
		skipK := true
		tolK := 1e-17
		tolu := 1e-12
		tols := 1e-14
		tests.CompareResults(tst, "data/spo751.sim", "cmp/spo751.cmp", "", tolK, tolu, tols, skipK, verb, nil)
	}
}

func Test_spo751re(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("spo751re. cylin exp. Richardson extrapolation")

	// run simulation
	main := fem.NewMain("data/spo751re.sim", "", true, false, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
	}

	// TODO: add check here
}
