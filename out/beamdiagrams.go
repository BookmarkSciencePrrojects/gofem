// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// DefineBeams define aliases for beams; e.g. "beam0", "beam1", etc.
func DefineBeams() {
	for _, beam := range Beams {
		alias := io.Sf("beam%d", beam.Id())
		Define(alias, P{{beam.Id(), -1}})
	}
}

// BeamDiagMoment plots beam's bending moment diagram
//  Input
//   alias -- use "" for all beams
//   idxI  -- index in TimeInds slice corresponding to selected output time; use -1 for the last item.
//           If alias defines a single point, the whole time series is returned and idxI is ignored.
//   withtext  -- show bending moment values
//   numfmt    -- number format for values. use "" for default
//   tolM      -- tolerance to clip absolute values of M
//   coef      -- coefficient to scale max(dimension) divided by max(Y); e.g. 0.1
//   sf        -- scaling factor. use 0 for automatic computation
func BeamDiagMoment(alias string, idxI int, withtext bool, numfmt string, tolM, coef, sf float64) {

	// all beams
	if alias == "" {

		// bending moments
		allM22 := make([][]float64, len(Beams))
		for i, beam := range Beams {
			alias = io.Sf("beam%d", beam.Id())
			allM22[i] = GetRes("M22", alias, idxI)
		}

		// scaling factor
		if sf < 1e-8 {
			maxAbsM := la.MatLargest(allM22, 1)
			dist := utl.Max(Dom.Msh.Xmax-Dom.Msh.Xmin, Dom.Msh.Ymax-Dom.Msh.Ymin)
			sf = 1.0
			if maxAbsM > 1e-7 {
				sf = coef * dist / maxAbsM
			}
		}

		// draw
		for i, beam := range Beams {
			beam.PlotDiagMoment(allM22[i], withtext, numfmt, tolM, sf)
		}
	}
}
