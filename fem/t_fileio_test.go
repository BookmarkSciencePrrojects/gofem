// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func Test_fileio01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("fileio01. File I/O operations")

	// start
	main := NewMain("data/bh16.sim", "", true, false, false, false, chk.Verbose, 0)

	// domain A
	domsA := NewDomains(main.Sim, main.DynCfs, 0, 1, false, false)
	if len(domsA) == 0 {
		tst.Errorf("NewDomains failed\n")
		return
	}
	domA := domsA[0]
	err := domA.SetStage(0)
	if err != nil {
		tst.Errorf("SetStage failed\n%v", err)
		return
	}
	for i, _ := range domA.Sol.Y {
		domA.Sol.Y[i] = float64(i)
	}
	io.Pforan("domA.Sol.Y = %v\n", domA.Sol.Y)

	// write file
	tidx := 123
	err = domA.SaveSol(tidx, true)
	if err != nil {
		tst.Errorf("SaveSol failed:\n%v", err)
		return
	}

	// domain B
	domsB := NewDomains(main.Sim, main.DynCfs, 0, 1, false, false)
	if len(domsB) == 0 {
		tst.Errorf("NewDomains failed\n")
		return
	}
	domB := domsB[0]
	err = domB.SetStage(0)
	if err != nil {
		tst.Errorf("SetStage failed\n%v", err)
		return
	}
	io.Pfpink("domB.Sol.Y (before) = %v\n", domB.Sol.Y)

	// read file
	err = domB.ReadSol(main.Sim.DirOut, main.Sim.Key, main.Sim.EncType, tidx)
	if err != nil {
		tst.Errorf("ReadSol failed:\n%v", err)
		return
	}
	io.Pfgreen("domB.Sol.Y (after) = %v\n", domB.Sol.Y)

	// check
	chk.Array(tst, "Y", 1e-17, domA.Sol.Y, domB.Sol.Y)
	chk.Array(tst, "dy/dt", 1e-17, domA.Sol.Dydt, domB.Sol.Dydt)
	chk.Array(tst, "d²y/dt²", 1e-17, domA.Sol.D2ydt2, domB.Sol.D2ydt2)
}
