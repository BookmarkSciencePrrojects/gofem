// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package diffusion

import (
	"testing"

	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/utl"
)

func Test_diffusion01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("diffusion01")

	// load sim => mesh => edat => cell
	sim := inp.ReadSim("data/diffusion.sim", "", true, true, 0)
	msh := sim.Regions[0].Msh
	edat := sim.Regions[0].ElemsData[0]
	cell := msh.Cells[0]

	// check info
	infofcn := ele.GetInfoFunc("diffusion")
	info := infofcn(sim, cell, edat)
	chk.IntAssert(len(info.Dofs), 9)
	for _, dof := range info.Dofs {
		chk.Strings(tst, "dofs", dof, []string{"u"})
	}
	chk.Strings(tst, "t1vars", info.T1vars, []string{"u"})

	// check element
	allocator := ele.GetAllocator("diffusion")
	e := allocator(sim, cell, edat, ele.BuildCoordsMatrix(cell, msh)).(*Diffusion)
	e.SetEqs([][]int{
		{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8},
	}, nil)
	chk.Ints(tst, "Umap", e.Umap, utl.IntRange(9))
}
