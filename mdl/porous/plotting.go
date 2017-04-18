// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package porous

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

// PlotLrm plots retention. If fname!="", figure is saved; otherwise it is not
//  Output:
//   X -- pc values
//   Y -- sl values
func PlotLrm(o *Model, dirout, fname string, pcmax float64, np int,
	returnTo0, withText, deriv bool,
	argsDry, argsWet, argsTxtMin, argsTxtMax, txtFmt string) (X, Y []float64) {

	// generate path
	if returnTo0 {
		X = utl.LinSpace(0, pcmax, np)
		X = append(X, utl.LinSpace(pcmax, 0, np)...)
	} else {
		X = utl.LinSpace(0, pcmax, np)
	}
	Y = make([]float64, len(X))
	var Z []float64
	if deriv {
		Z = make([]float64, len(X)) // dsl/dpc
	}

	// allocate auxiliary state variables
	ρL, ρG, pl, pg := o.Liq.R0, o.Gas.R0, 0.0, 0.0
	sta, err := o.NewState(ρL, ρG, pl, pg)
	if err != nil {
		chk.Panic("cannot create new state:\n%v", err)
	}

	// compute LRM
	Y[0] = o.Lrm.SlMax()
	for i := 1; i < len(X); i++ {
		Δpc := X[i] - X[i-1]
		pl := -X[i]
		err = o.Update(sta, -Δpc, 0, pl, 0)
		if err != nil {
			chk.Panic("cannot update state:\n%v", err)
		}
		Y[i] = sta.A_sl
		if deriv {
			Z[i], err = o.Ccb(sta, X[i])
			if err != nil {
				chk.Panic("cannot compute Ccb:\n%v", err)
			}
		}
	}

	// plot LRM
	if deriv {
		plt.Subplot(2, 1, 1)
	}
	if argsDry == "" {
		argsDry = "'b-', clip_on=0, label='drying'"
	}
	if argsWet == "" {
		argsWet = "'b-', clip_on=0, label='wetting', color='#0397dc'"
	}
	//if returnTo0 {
	//plt.Plot(X[:np], Y[:np], argsDry)
	//plt.Plot(X[np:], Y[np:], argsWet)
	//} else {
	//plt.Plot(X, Y, argsDry)
	//}

	// add text
	if withText {
		//l := np - 1
		if argsTxtMin == "" {
			argsTxtMin = "ha='left',  color='red', size=8"
		}
		if argsTxtMax == "" {
			argsTxtMax = "ha='right', color='red', size=8"
		}
		if txtFmt == "" {
			txtFmt = "%g"
		}
		//plt.Text(X[0], Y[0], io.Sf("(%g, %g)", X[0], Y[0]), argsTxtMin)
		//plt.Text(X[l], Y[l], io.Sf("(%g, "+txtFmt+")", X[l], Y[l]), argsTxtMax)
	}
	//plt.Gll("$p_c$", "$s_{\\ell}$", "")

	// plot deriatives
	if deriv {
		plt.Subplot(2, 1, 2)
		//if returnTo0 {
		//plt.Plot(X[:np], Z[:np], "'b-', clip_on=0, color='#0397dc'")
		//plt.Plot(X[np:], Z[np:], "'b-', clip_on=0")
		//} else {
		//plt.Plot(X, Y, "'b-', clip_on=0")
		//}
		//if withText {
		//l := np - 1
		//plt.Text(X[0], Z[0], io.Sf("(%g, %g)", X[0], Z[0]), "ha='left',  color='red', size=8")
		//plt.Text(X[l], Z[l], io.Sf("(%g, %g)", X[l], Z[l]), "ha='right', color='red', size=8")
		//}
		//plt.Gll("$p_c$", "$\\bar{C}_c=\\mathrm{d}s_{\\ell}/\\mathrm{d}p_c$", "")
	}

	// save figure
	if fname != "" {
		plt.SaveD(dirout, fname)
	}
	return
}
