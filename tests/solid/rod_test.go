// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"testing"

	"github.com/cpmech/gofem/ele/solid"
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/tests"
	"github.com/cpmech/gosl/chk"
)

func Test_bridge01a(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("bridge01a. simple bridge section")

	// fem
	main := fem.NewMain("data/bridge01.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 1e-11
	tolu := 1e-15
	tols := 1e-9
	tests.CompareResults(tst, "data/bridge01.sim", "cmp/bridge01.cmp", "", tolK, tolu, tols, skipK, chk.Verbose, nil)
}

func Test_bridge01b(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("bridge01. simple bridge section. ElastRod")

	// fem
	main := fem.NewMain("data/bridge01erod.sim", "", true, true, false, false, chk.Verbose, 0)

	// set stage
	err := main.SetStage(0)
	if err != nil {
		tst.Errorf("SetStage failed:\n%v", err)
		return
	}

	// recompute matrices
	for _, elem := range main.Domains[0].Elems {
		e := elem.(*solid.ElastRod)
		e.Recompute(true)
	}

	// run
	err = main.SolveOneStage(0, true)
	if err != nil {
		tst.Error("SolveOneStage failed:\n%v", err)
		return
	}

	// check
	skipK := false
	tolK := 1e-11
	tolu := 1e-15
	tols := 1e-9
	tests.CompareResults(tst, "data/bridge01erod.sim", "cmp/bridge01.cmp", "", tolK, tolu, tols, skipK, chk.Verbose, nil)
}
