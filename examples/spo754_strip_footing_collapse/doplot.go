// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"testing"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/plt"
)

func main() {

	// options
	verbose := false

	// start analysis process
	out.Start("spo754.sim", 0, 0)

	// all nodes
	out.Define("all nodes", out.AllNodes())
	out.Define("A", out.At{0, 5})

	// load results
	out.LoadResults(nil)

	// check
	skipK := true
	tolK := 1e-17
	tolu := 1e-11
	tols := 1e-04

	var tst testing.T
	fem.TestingCompareResultsU(&tst, "spo754.sim", "spo754.cmp", "", tolK, tolu, tols, skipK, verbose)

	// save
	plt.SetForPng(0.8, 400, 200)
	out.Plot("uy", "t", "A", plt.Fmt{C: "r", M: "o"}, -1)
	out.Draw("/tmp", "spo754.png", -1, -1, false, nil)
}
