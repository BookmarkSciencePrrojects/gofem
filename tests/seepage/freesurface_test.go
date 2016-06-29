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
)

func Test_frees01a(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("frees01a. Freesurface liquid-only. Check DOFs and BCs")

	// start simulation
	main := fem.NewMain("data/frees01.sim", "", true, false, false, false, chk.Verbose, 0)

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
	chk.IntAssert(len(dom.Nodes), 62)
	chk.IntAssert(len(dom.Elems), 15)

	// vertices with "fl"
	seepverts := map[int]bool{3: true, 45: true, 7: true, 49: true, 11: true, 53: true, 15: true, 57: true, 19: true, 61: true, 23: true}

	// check dofs
	var seepeqs []int
	for _, nod := range dom.Nodes {
		if seepverts[nod.Vert.Id] {
			chk.IntAssert(len(nod.Dofs), 2)
			seepeqs = append(seepeqs, nod.Dofs[1].Eq)
		} else {
			chk.IntAssert(len(nod.Dofs), 1)
		}
	}
	sort.Ints(seepeqs)
	io.Pforan("seepeqs = %v\n", seepeqs)
	chk.Ints(tst, "seepeqs", seepeqs, []int{14, 16, 19, 30, 32, 43, 45, 56, 58, 69, 71})

	// check Fmap
	e2 := dom.Elems[2].(*seepage.Liquid)
	chk.Ints(tst, "e2.Fmap", e2.Fmap, []int{14, 16, 19})
	e5 := dom.Elems[5].(*seepage.Liquid)
	chk.Ints(tst, "e5.Fmap", e5.Fmap, []int{16, 30, 32})
	e8 := dom.Elems[8].(*seepage.Liquid)
	chk.Ints(tst, "e8.Fmap", e8.Fmap, []int{30, 43, 45})
	e11 := dom.Elems[11].(*seepage.Liquid)
	chk.Ints(tst, "e11.Fmap", e11.Fmap, []int{43, 56, 58})
	e14 := dom.Elems[14].(*seepage.Liquid)
	chk.Ints(tst, "e14.Fmap", e14.Fmap, []int{56, 69, 71})
}

func Test_frees01b(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("frees01b. Freesurface liquid-only. Run and check Kb")

	// start simulation
	main := fem.NewMain("data/frees01.sim", "", true, true, false, false, chk.Verbose, 0)

	// for debugging Kb
	if true {
		tests.Liquid(main, &tests.Kb{
			Tst: tst, Eid: 14, Tol: 1e-5, Verb: chk.Verbose,
			Ni: 1, Nj: 1, ItMin: 1, ItMax: -1, Tmin: 200, Tmax: 200,
		})
	}

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
	}
}
