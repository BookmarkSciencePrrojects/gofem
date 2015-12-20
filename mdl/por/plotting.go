// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package por

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func PlotSimple(o *Model, dirout, fname string, pcmax float64, np int, returnTo0, withText, deriv bool) {
	var X []float64 // pc
	var Y []float64 // sl
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
	ρL, ρG, pl, pg := o.RhoL0, o.RhoG0, 0.0, 0.0
	sta, err := o.NewState(ρL, ρG, pl, pg)
	if err != nil {
		chk.Panic("cannot create new state:\n%v", err)
	}
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
	if deriv {
		plt.Subplot(2, 1, 1)
	}
	if returnTo0 {
		plt.Plot(X[:np], Y[:np], "'b-', clip_on=0, color='#0397dc'")
		plt.Plot(X[np:], Y[np:], "'b-', clip_on=0")
	} else {
		plt.Plot(X, Y, "'b-', clip_on=0")
	}
	if withText {
		l := np - 1
		plt.Text(X[0], Y[0], io.Sf("(%g, %g)", X[0], Y[0]), "ha='left',  color='red', size=8")
		plt.Text(X[l], Y[l], io.Sf("(%g, %g)", X[l], Y[l]), "ha='right', color='red', size=8")
	}
	plt.Gll("$p_c$", "$s_{\\ell}$", "")
	if deriv {
		plt.Subplot(2, 1, 2)
		if returnTo0 {
			plt.Plot(X[:np], Z[:np], "'b-', clip_on=0, color='#0397dc'")
			plt.Plot(X[np:], Z[np:], "'b-', clip_on=0")
		} else {
			plt.Plot(X, Y, "'b-', clip_on=0")
		}
		if withText {
			l := np - 1
			plt.Text(X[0], Z[0], io.Sf("(%g, %g)", X[0], Z[0]), "ha='left',  color='red', size=8")
			plt.Text(X[l], Z[l], io.Sf("(%g, %g)", X[l], Z[l]), "ha='right', color='red', size=8")
		}
		plt.Gll("$p_c$", "$\\bar{C}_c=\\mathrm{d}s_{\\ell}/\\mathrm{d}p_c$", "")
	}
	plt.SaveD(dirout, fname)
}
