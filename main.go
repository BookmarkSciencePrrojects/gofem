// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"flag"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/mpi"
	"github.com/cpmech/gosl/utl"
)

func main() {

	// in put data
	erasefiles := true
	verbose := true

	// catch errors
	defer func() {
		if err := recover(); err != nil {
			if mpi.Rank() == 0 {
				chk.Verbose = true
				for i := 8; i > 3; i-- {
					chk.CallerInfo(i)
				}

				// print log file
				if verbose {
					buf, err := io.ReadFile(inp.LogFile.Name())
					if err != nil {
						io.Pfred("cannot read log file:%v\n", err)
					}
					io.Pfyel("\n%v\n", string(buf))
				}

				io.PfRed("ERROR: %v\n", err)
			}
		}
		// make sure to flush log
		defer fem.End()
		mpi.Stop(false)
	}()
	mpi.Start(false)

	// message
	if mpi.Rank() == 0 {
		io.PfWhite("\nGofem v3 -- Go Finite Element Method\n\n")
		io.Pf("Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.\n")
		io.Pf("Use of this source code is governed by a BSD-style\n")
		io.Pf("license that can be found in the LICENSE file.\n\n")
	}

	// simulation filenamepath
	flag.Parse()
	var fnamepath string
	if len(flag.Args()) > 0 {
		fnamepath = flag.Arg(0)
	} else {
		chk.Panic("Please, provide a filename. Ex.: cylinder.sim")
	}

	// check extension
	if io.FnExt(fnamepath) == "" {
		fnamepath += ".sim"
	}

	// other options
	if len(flag.Args()) > 1 {
		erasefiles = io.Atob(flag.Arg(1))
	}
	if len(flag.Args()) > 2 {
		verbose = io.Atob(flag.Arg(2))
	}

	// profiling?
	defer utl.DoProf(false)()

	// start global variables and log
	if !fem.Start(fnamepath, erasefiles, verbose) {
		chk.Panic("Start failed\n")
		return
	}

	// run simulation
	if !fem.Run() {
		chk.Panic("Run failed\n")
	}
}
