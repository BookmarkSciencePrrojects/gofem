// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"testing"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/plt"
)

func Test_beamdiag01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("beamdiag01")

	// run only if verbose==true
	if chk.Verbose {

		// start simulation
		main := fem.NewMain("../fem/data/beam03.sim", "", true, true, false, false, chk.Verbose, 0)

		// set stage
		stgidx := 2
		err := main.SetStage(stgidx)
		if err != nil {
			tst.Errorf("SetStage failed:\n%v", err)
			return
		}

		// run one stage
		err = main.SolveOneStage(stgidx, true)
		if err != nil {
			tst.Errorf("SolveOneStage failed:\n%v", err)
			return
		}

		// start post-processing
		Start("../fem/data/beam03.sim", stgidx, 0)

		// define points
		DefineBeams()

		// load results
		LoadResults(nil)

		// plot bending moment diagram
		withtext, numfmt, tol, coef, sf, onlyLin := true, "", 1e-10, 0.2, 0.0, false
		plt.SetForPng(1, 600, 150)
		Dom.Msh.Draw2d(onlyLin, false, nil, 3)
		BeamDiagMoment("", -1, withtext, numfmt, tol, coef, sf)
		plt.SaveD("/tmp/gofem", "test_beamdiag02.png")
	}
}
