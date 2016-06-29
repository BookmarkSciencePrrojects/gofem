// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package porous

import (
	"testing"

	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/chk"
)

func Test_upp01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("upp01")

	// load sim => mesh => edat => cell
	sim := inp.ReadSim("data/solid-liquid-gas.sim", "", true, 0)
	msh := sim.Regions[0].Msh
	edat := sim.Regions[0].ElemsData[0]
	cell := msh.Cells[0]

	// check info
	infofcn := ele.GetInfoFunc("solid-liquid-gas")
	info := infofcn(sim, cell, edat)
	chk.IntAssert(len(info.Dofs), 9)
	for m, dof := range info.Dofs {
		if m < 4 {
			chk.Strings(tst, "upp dofs", dof, []string{"ux", "uy", "pl", "pg"})
		} else {
			chk.Strings(tst, "u dofs", dof, []string{"ux", "uy"})
		}
	}
	chk.Strings(tst, "t1vars", info.T1vars, []string{"pl", "pg"})
	chk.Strings(tst, "t2vars", info.T2vars, []string{"ux", "uy"})

	// check element
	allocator := ele.GetAllocator("solid-liquid-gas")
	e := allocator(sim, cell, edat, ele.BuildCoordsMatrix(cell, msh)).(*SolidLiquidGas)
	e.SetEqs([][]int{
		{0, 1, 2, 3},
		{4, 5, 6, 7},
		{8, 9, 10, 11},
		{12, 13, 14, 15},
		{16, 17},
		{18, 19},
		{20, 21},
		{22, 23},
		{24, 25},
	}, nil)
	chk.Ints(tst, "Umap", e.U.Umap, []int{0, 1, 4, 5, 8, 9, 12, 13, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25})
}
