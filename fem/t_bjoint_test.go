// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func Test_bjoint01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("bjoint01. beam joint compatible. pull-out")

	// start simulation
	analysis := NewFEM("data/bjointcomp2d01.sim", "", true, false, false, false, chk.Verbose, 0)

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
	chk.IntAssert(len(dom.Nodes), 12)
	chk.IntAssert(len(dom.Elems), 7)
}
