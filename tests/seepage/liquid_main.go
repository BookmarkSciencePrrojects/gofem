// Copyright 2016 The Gofem Authors. All rights reserved.
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
	main := fem.NewMain("data/p01.sim", "", true, true, false, true, true, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Error("Run failed\n")
		return
	}
}
