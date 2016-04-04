// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"
	"sort"
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func Test_diffu01a(tst *testing.T) {

	/* Diffusion equation
	 *
	 *      Nodes / Tags                 Equations
	 *
	 *     8 o----o----o 9 (-5)      22 o----o----o 21
	 *       |   14    |                |   24    |
	 *       |         |                |         |
	 *    21 o    o    o 22 (-6)     25 o    o    o 23
	 *       |   26    |                |   26    |
	 *       |         |                |         |
	 *     6 o----o----o 7 (-4)      16 o----o----o 15
	 *       |   13    |                |   18    |
	 *       |         |                |         |
	 *    19 o    o    o 20 (-6)     19 o    o    o 17
	 *       |   25    |                |   20    |
	 *       |         |                |         |
	 *     4 o----o----o 5 (-3)      10 o----o----o 9
	 *       |   12    |                |   12    |
	 *       |         |                |         |
	 *    17 o    o    o 18 (-6)     13 o    o    o 11
	 *       |   24    |                |   14    |
	 *       |         |                |         |
	 *     2 o----o----o 3 (-2)       3 o----o----o 2
	 *       |   11    |                |    6    |
	 *       |         |                |         |
	 *    15 o    o    o 16 (-6)      7 o    o    o 5
	 *       |   23    |                |    8    |
	 *       |         |                |         |
	 *     0 o----o----o 1 (-1)       0 o----o----o 1
	 *           10                          4
	 */

	//verbose()
	chk.PrintTitle("diffu01a. Poisson equation 01. Check DOFs")

	// start simulation
	analysis := NewFEM("data/diffu01.sim", "", true, false, false, false, chk.Verbose, 0)

	// set stage
	err := analysis.SetStage(0)
	if err != nil {
		tst.Errorf("SetStage failed:\n%v", err)
		return
	}

	// initialise solution vectors
	err = analysis.ZeroStage(0, true)
	if err != nil {
		tst.Errorf("ZeroStage failed:\n%v", err)
		return
	}

	// domain
	dom := analysis.Domains[0]

	// nodes and elements
	chk.IntAssert(len(dom.Nodes), 27)
	chk.IntAssert(len(dom.Elems), 4)

	// check dofs
	for _, nod := range dom.Nodes {
		chk.IntAssert(len(nod.Dofs), 1)
		chk.StrAssert(nod.Dofs[0].Key, "u")
	}

	// check equations
	nids, eqs := get_nids_eqs(dom)
	chk.Ints(tst, "nids", nids, []int{
		0, 1, 3, 2, 10, 16, 11, 15, 23,
		5, 4, 18, 12, 17, 24,
		7, 6, 20, 13, 19, 25,
		9, 8, 22, 14, 21, 26,
	})
	chk.Ints(tst, "eqs", eqs, utl.IntRange(27))

	// check Umap
	Umaps := [][]int{
		{0, 1, 2, 3, 4, 5, 6, 7, 8},
		{3, 2, 9, 10, 6, 11, 12, 13, 14},
		{10, 9, 15, 16, 12, 17, 18, 19, 20},
		{16, 15, 21, 22, 18, 23, 24, 25, 26},
	}
	for i, ele := range dom.Elems {
		e := ele.(*ElemDiffu)
		io.Pforan("e%d.Umap = %v\n", e.Id(), e.Umap)
		chk.Ints(tst, "Umap", e.Umap, Umaps[i])
	}

	// constraints
	chk.IntAssert(len(dom.EssenBcs.Bcs), 6)
	var ct_u_eqs []int // equations with u prescribed [sorted]
	for _, c := range dom.EssenBcs.Bcs {
		chk.IntAssert(len(c.Eqs), 1)
		eq := c.Eqs[0]
		io.Pforan("key=%v eq=%v\n", c.Key, eq)
		switch c.Key {
		case "u":
			ct_u_eqs = append(ct_u_eqs, eq)
		default:
			tst.Errorf("key %s is incorrect", c.Key)
		}
	}
	sort.Ints(ct_u_eqs)
	chk.Ints(tst, "equations with u prescribed", ct_u_eqs, []int{0, 1, 4, 21, 22, 24})
}

func Test_diffu01b(tst *testing.T) {

	//verbose()
	chk.PrintTitle("diffu01b. Poisson equation 01")
	doplot := true
	if !chk.Verbose {
		doplot = false
	}

	// run simulation
	analysis := NewFEM("data/diffu01.sim", "", true, false, false, false, chk.Verbose, 0)

	// for debugging Kb
	if true {
		p_DebugKb(analysis, &testKb{
			tst: tst, eid: 3, tol: 1e-6, verb: chk.Verbose,
			ni: -1, nj: -1, itmin: -1, itmax: -1,
		})
	}

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// analytical solution
	L := 10.0
	C := L * L / 6.0
	uanaF := func(x []float64) float64 {
		return -math.Pow(x[1], 3.0)/6.0 + C*x[1]
	}

	// check
	dom := analysis.Domains[0]
	var Y, Unum, Uana []float64
	for _, nod := range dom.Nodes {
		x := nod.Vert.C
		eq := nod.GetEq("u")
		unum := dom.Sol.Y[eq]
		uana := uanaF(x)
		chk.AnaNum(tst, io.Sf("u(%6.3f)", x), 1e-12, unum, uana, chk.Verbose)
		if doplot {
			if math.Abs(x[0]) < 1e-10 {
				Y = append(Y, x[1])
				Unum = append(Unum, unum)
				Uana = append(Uana, uana)
			}
		}
	}
	if doplot {
		_, Y, Unum, Uana, err = utl.SortQuadruples(nil, Y, Unum, Uana, "x")
		if err != nil {
			tst.Errorf("cannot sort arrays: %v\n", err)
			return
		}
		plt.Plot(Y, Unum, "'ks', label='numerical'")
		plt.Plot(Y, Uana, "'r-', label='analytical'")
		plt.Gll("y", "u", "")
		plt.SaveD("/tmp/gofem", "test_diffu01b.png")
	}
}
