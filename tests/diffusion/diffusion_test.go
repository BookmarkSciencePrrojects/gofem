// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"math"
	"sort"
	"testing"

	"github.com/cpmech/gofem/ele/diffusion"
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/tests"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func Test_diffu01a(tst *testing.T) {

	/* Diffusion equation
	 *
	 *        Nodes                Equations
	 *
	 *     0     1     2         3     6     2
	 *     o-----o-----o         o-----o-----o
	 *     |           |         |           |
	 *     |     4     |         |     8     |
	 *   3 o     o     o 5     7 o     o     o 5
	 *     |           |         |           |
	 *     |           |         |           |
	 *     o-----o-----o         o-----o-----o
	 *     6     7     8         0     4     1
	 */

	//tests.Verbose()
	chk.PrintTitle("diffu01a. Diffusion (Poisson) equation 01. Check DOFs")

	// start simulation
	main := fem.NewMain("data/diffu01.sim", "", true, false, false, false, chk.Verbose, 0)

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
	io.Pforan("dom.elems = %+#v\n", dom.Elems)

	// nodes and elements
	chk.IntAssert(len(dom.Nodes), 9)
	chk.IntAssert(len(dom.Elems), 1)

	// check dofs
	for _, nod := range dom.Nodes {
		chk.IntAssert(len(nod.Dofs), 1)
	}

	// check equations
	nids, eqs := tests.GetNidsEqs(dom)
	chk.Ints(tst, "nids", nids, []int{6, 8, 2, 0, 7, 5, 1, 3, 4})
	chk.Ints(tst, "eqs", eqs, []int{0, 1, 2, 3, 4, 5, 6, 7, 8})

	// check solution arrays
	ny := 9
	nλ := 5
	nyb := ny + nλ
	chk.IntAssert(len(dom.Sol.Y), ny)
	chk.IntAssert(len(dom.Sol.Dydt), 0)
	chk.IntAssert(len(dom.Sol.L), nλ)
	chk.IntAssert(len(dom.Sol.ΔY), ny)

	// check linear solver arrays
	chk.IntAssert(len(dom.Fb), nyb)
	chk.IntAssert(len(dom.Wb), nyb)

	// check Tmap
	Umaps := [][]int{ // for each element
		{0, 1, 2, 3, 4, 5, 6, 7, 8}, // element # 0
	}
	for i, element := range dom.Elems {
		e := element.(*diffusion.Diffusion)
		io.Pforan("e%d.Umap = %v\n", e.Id(), e.Umap)
		chk.Ints(tst, "Umap", e.Umap, Umaps[i])
	}

	// constraints
	chk.IntAssert(len(dom.EssenBcs.Bcs), nλ)
	var ct_u_eqs []int // constrained T equations [sorted]
	for _, c := range dom.EssenBcs.Bcs {
		io.Pforan("c = %+#v\n", c)
		chk.IntAssert(len(c.Eqs), 1)
		chk.IntAssert(len(c.ValsA), 1)
		chk.String(tst, c.Key, "u")
		chk.Scalar(tst, "A", 1e-15, c.ValsA[0], 1.0)
		eq := c.Eqs[0]
		ct_u_eqs = append(ct_u_eqs, eq)
	}
	sort.Ints(ct_u_eqs)
	chk.Ints(tst, "constrained u equations", ct_u_eqs, []int{0, 2, 3, 6, 7})
}

func Test_diffu02a(tst *testing.T) {

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

	//tests.Verbose()
	chk.PrintTitle("diffu02a. Diffusion (Poisson) equation 02. Check DOFs")

	// start simulation
	main := fem.NewMain("data/diffu02.sim", "", true, false, false, false, chk.Verbose, 0)

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
		chk.IntAssert(len(nod.Dofs), 1)
		chk.StrAssert(nod.Dofs[0].Key, "u")
	}

	// check equations
	nids, eqs := tests.GetNidsEqs(dom)
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
		e := ele.(*diffusion.Diffusion)
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

func Test_diffu02b(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("diffu02b. Diffusion (Poisson) equation 02. Run")
	doplot := true
	if !chk.Verbose {
		doplot = false
	}

	// run simulation
	main := fem.NewMain("data/diffu02.sim", "", true, false, false, false, chk.Verbose, 0)

	// for debugging Kb
	if true {
		tests.Liquid(main, &tests.Kb{
			Tst: tst, Eid: 3, Tol: 1e-6, Verb: chk.Verbose,
			Ni: -1, Nj: -1, ItMin: -1, ItMax: -1,
		})
	}

	// run simulation
	err := main.Run()
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
	dom := main.Domains[0]
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

func test_diffu03(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("diffu03. Diffusion (Poisson) equation 03. Refinement")

	// start simulation
	main := fem.NewMain("data/diffu03.sim", "", true, false, false, false, chk.Verbose, 0)

	// refine mesh
	dom := main.Domains[0]
	msh := dom.Msh
	io.Pforan("msh = %v\n", msh)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// analytical solution
	w, H, Ttop, Tside := 1.0, 3.0, 60.0, 20.0
	pi := math.Pi
	npts := 100
	calc_u := func(x []float64) float64 {
		sum := 0.0
		for i := 1; i < npts; i++ {
			n := float64(i)
			sum += ((math.Pow(-1, n+1) + 1) / n) * math.Sin(n*pi*x[0]/w) * math.Sinh(n*pi*x[1]/w) / math.Sinh(n*pi*H/w)
		}
		return (Ttop-Tside)*(2.0/pi)*sum + Tside
	}

	// check
	for _, nod := range dom.Nodes {
		x := nod.Vert.C
		eq := nod.GetEq("u")
		unum := dom.Sol.Y[eq]
		uana := calc_u(x)
		io.Pf("%2d : u(%6.3f) = %23g   %23g\n", nod.Vert.Id, x, unum, uana)
		//chk.AnaNum(tst, io.Sf("u(%6.3f)", x), 1e-12, unum, uana, chk.Verbose)
	}
}
