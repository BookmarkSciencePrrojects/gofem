// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package conduct

import (
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func Plot(o Model, dirout, fname string, np int, gas, withText, deriv bool) {
	X := utl.LinSpace(0, 1, np)
	Y := make([]float64, np)
	var Z []float64
	if deriv {
		Z = make([]float64, np)
	}
	for i := 0; i < np; i++ {
		if gas {
			Y[i] = o.Kgr(X[i])
		} else {
			Y[i] = o.Klr(X[i])
		}
		if deriv {
			if gas {
				Z[i] = o.DkgrDsg(X[i])
			} else {
				Z[i] = o.DklrDsl(X[i])
			}
		}
	}
	if deriv {
		plt.Subplot(2, 1, 1)
	}
	plt.Plot(X, Y, "'b-', clip_on=0")
	if withText {
		l := np - 1
		plt.Text(X[0], Y[0], io.Sf("(%g, %g)", X[0], Y[0]), "ha='left',  color='red', size=8")
		plt.Text(X[l], Y[l], io.Sf("(%g, %g)", X[l], Y[l]), "ha='right', color='red', size=8")
	}
	key := "\\ell"
	if gas {
		key = "g"
	}
	plt.Gll("$s_{"+key+"}$", "$k_{"+key+"}^r$", "")
	if deriv {
		plt.Subplot(2, 1, 2)
		plt.Plot(X, Z, "'b-', clip_on=0")
		if withText {
			l := np - 1
			plt.Text(X[0], Z[0], io.Sf("(%g, %g)", X[0], Z[0]), "ha='left',  color='red', size=8")
			plt.Text(X[l], Z[l], io.Sf("(%g, %g)", X[l], Z[l]), "ha='right', color='red', size=8")
		}
		plt.Gll("$s_{"+key+"}$", "$\\mathrm{d}{k_{"+key+"}^r}/\\mathrm{d}{s_{"+key+"}}$", "")
	}
	plt.SaveD(dirout, fname)
}
