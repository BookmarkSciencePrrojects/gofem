// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package seepage

import (
	"testing"

	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/utl"
)

func Test_liquid01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("liquid01")

	// load sim => mesh => edat => cell
	sim := inp.ReadSim("data/liquid.sim", "", true, 0)
	msh := sim.Regions[0].Msh
	edat := sim.Regions[0].ElemsData[0]
	cell := msh.Cells[0]

	// check info
	infofcn := ele.GetInfoFunc("liquid")
	info := infofcn(sim, cell, edat)
	chk.IntAssert(len(info.Dofs), 9)
	for _, dof := range info.Dofs {
		chk.Strings(tst, "pl dofs", dof, []string{"pl"})
	}
	chk.Strings(tst, "t1vars", info.T1vars, []string{"pl"})

	// check element
	allocator := ele.GetAllocator("liquid")
	e := allocator(sim, cell, edat, ele.BuildCoordsMatrix(cell, msh)).(*Liquid)
	e.SetEqs([][]int{
		{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8},
	}, nil)
	chk.Ints(tst, "Pmap", e.Pmap, utl.IntRange(9))
}
