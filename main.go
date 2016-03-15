// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/mpi"
	"github.com/cpmech/gosl/utl"
)

func main() {

	// catch errors
	defer func() {
		if err := recover(); err != nil {
			if mpi.Rank() == 0 {
				io.PfRed("\nERROR: %v", err)
				io.Pf("See location of error below:\n")
				chk.Verbose = true
				for i := 5; i > 3; i-- {
					chk.CallerInfo(i)
				}
			}
		}
		mpi.Stop(false)
	}()
	mpi.Start(false)

	// read input parameters
	fnamepath, _ := io.ArgToFilename(0, "", ".sim", true)
	verbose := io.ArgToBool(1, true)
	erasePrev := io.ArgToBool(2, true)
	saveSummary := io.ArgToBool(3, true)
	allowParallel := io.ArgToBool(4, true)
	doprof := io.ArgToInt(5, 0)

	// message
	if mpi.Rank() == 0 && verbose {
		io.PfWhite("\nGofem Version 3.1 -- Go Finite Element Method\n")
		io.Pf("Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.\n")
		io.Pf("Use of this source code is governed by a BSD-style\n")
		io.Pf("license that can be found in the LICENSE file.\n")

		io.Pf("\n%v\n", io.ArgsTable("INPUT ARGUMENTS",
			"filename path", "fnamepath", fnamepath,
			"show messages", "verbose", verbose,
			"erase previous results", "erasePrev", erasePrev,
			"save summary", "saveSummary", saveSummary,
			"allow parallel run", "allowParallel", allowParallel,
			"profiling: 0=none 1=CPU 2=MEM", "doprof", doprof,
		))
	}

	// profiling?
	if doprof > 0 {
		defer utl.DoProf(false, doprof)()
	}

	// analysis data
	alias := ""
	readSummary := false
	analysis := fem.NewFEM(fnamepath, alias, erasePrev, saveSummary, readSummary, allowParallel, verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		chk.Panic("Run failed:\n%v", err)
	}
}
