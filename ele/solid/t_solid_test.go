// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package solid

import (
	"testing"

	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/utl"
)

func Test_solid01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("solid01")

	// load sim => mesh => edat => cell
	sim := inp.ReadSim("data/solid.sim", "", true, 0)
	msh := sim.Regions[0].Msh
	edat := sim.Regions[0].ElemsData[0]
	cell := msh.Cells[0]

	// check info
	infofcn := ele.GetInfoFunc("solid")
	info := infofcn(sim, cell, edat)
	chk.IntAssert(len(info.Dofs), 9)
	for _, dof := range info.Dofs {
		chk.Strings(tst, "u dofs", dof, []string{"ux", "uy"})
	}
	chk.Strings(tst, "t2vars", info.T2vars, []string{"ux", "uy"})

	// check element
	allocator := ele.GetAllocator("solid")
	e := allocator(sim, cell, edat, ele.BuildCoordsMatrix(cell, msh)).(*Solid)
	e.SetEqs([][]int{
		{0, 1}, {2, 3}, {4, 5}, {6, 7}, {8, 9}, {10, 11}, {12, 13}, {14, 15}, {16, 17},
	}, nil)
	chk.Ints(tst, "Umap", e.Umap, utl.IntRange(18))
}
