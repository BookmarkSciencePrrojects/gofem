// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"testing"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/mpi"
)

func main() {

	// catch errors
	var tst testing.T
	defer func() {
		if mpi.Rank() == 0 {
			if err := recover(); err != nil {
				io.PfRed("ERROR: %v\n", err)
			}
			if tst.Failed() {
				io.PfRed("test failed\n")
			}
		}
		mpi.Stop(false)
	}()
	mpi.Start(false)

	// start global variables and log
	analysis := fem.NewFEM("data/spo751.sim", "", true, true, false, true, true, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Error("Run failed\n")
		return
	}

	// check
	skipK := true
	tolK := 1e-17
	tolu := 1e-12
	tols := 1e-14
	fem.TestingCompareResultsU(&tst, "data/spo751.sim", "cmp/spo751.cmp", "", tolK, tolu, tols, skipK, true)
}
