// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// PlotAllBendingMom2d plots all bending moments (2D beams)
//  Input:
//   dom       -- Domain
//   nstations -- number of stations
//   withtext  -- show bending moment values
//   numfmt    -- number format for values. use "" for default
//   tolM      -- tolerance to clip absolute values of M
//   coef      -- coefficient to scale max(dimension) divided by max(Y); e.g. 0.1
//   sf        -- scaling factor. use 0 for automatic computation
//  Output:
//   beams  -- all beam elements
//   allMrr -- all Mrr bending moments at all beams
func PlotAllBendingMom2d(dom *Domain, nstations int, withtext bool, numfmt string, tolM, coef, sf float64) (beams []*Beam, allMrr [][]float64) {

	// collect beams
	for _, elem := range dom.Elems {
		if beam, ok := elem.(*Beam); ok {
			beams = append(beams, beam)
		}
	}

	// compute bending moments
	allMrr = make([][]float64, len(beams))
	for i, beam := range beams {
		allMrr[i] = beam.CalcMoment2d(dom.Sol, 0, nstations)
	}

	// scaling factor
	if sf < 1e-8 {
		maxAbsM := la.MatLargest(allMrr, 1)
		dist := utl.Max(dom.Msh.Xmax-dom.Msh.Xmin, dom.Msh.Ymax-dom.Msh.Ymin)
		sf = 1.0
		if maxAbsM > 1e-7 {
			sf = coef * dist / maxAbsM
		}
	}

	// draw
	for i, beam := range beams {
		beam.PlotDiagMoment(allMrr[i], withtext, numfmt, tolM, sf)
	}
	return
}
