// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func main() {

	// filename
	filename, fnkey := io.ArgToFilename(0, "liqgas-3mcol-4el", ".sim", true)
	fourE := fnkey == "liqgas-3mcol-4el"

	// start analysis process
	out.Start(filename, 0, 0)

	// define entities
	out.Define("A(top) B C(mid) D E(bot)", out.N{-5, -4, -3, -2, -1})
	if fourE {
		out.Define("a(top) b c d(bot)", out.P{{3, 6}, {2, 3}, {1, 3}, {0, 0}})
	} else {
		out.Define("a(top) b c(mid) d e(bot)", out.P{{15, 6}, {12, 3}, {8, 0}, {4, 3}, {0, 0}})
	}
	out.Define("left-side", out.Along{{0, 0}, {0, 10}})

	// load results
	out.LoadResults(nil)

	// styles
	me := 10
	S := []plt.Fmt{
		plt.Fmt{C: "b", M: "*", Me: me},
		plt.Fmt{C: "g", M: "o", Me: me},
		plt.Fmt{C: "m", M: "x", Me: me},
		plt.Fmt{C: "orange", M: "+", Me: me},
		plt.Fmt{C: "r", M: "^", Me: me},
	}

	// pl
	out.Splot("liquid pressure along column")
	I, _ := utl.GetITout(out.Times, []float64{0, 500, 1000, 1500, 2000}, 0.1)
	for _, i := range I {
		t := out.Times[i]
		out.Plot("pl", "y", "left-side", plt.Fmt{L: io.Sf("t=%g", t)}, i)
	}

	// pg
	out.Splot("gas pressure along column")
	for _, i := range I {
		t := out.Times[i]
		out.Plot("pg", "y", "left-side", plt.Fmt{L: io.Sf("t=%g", t)}, i)
	}

	// pl
	out.Splot("liquid pressure versus time")
	for i, l := range []string{"A(top)", "B", "C(mid)", "D", "E(bot)"} {
		out.Plot("t", "pl", l, S[i], -1)
	}
	//out.Csplot.GllArgs = "leg_out=1, leg_ncol=5, leg_hlen=2"

	// pg
	out.Splot("gas pressure versus time")
	for i, l := range []string{"A(top)", "B", "C(mid)", "D", "E(bot)"} {
		out.Plot("t", "pg", l, S[i], -1)
	}

	// sl
	out.Splot("liquid saturation")
	ekeys := []string{"a(top)", "b", "c(mid)", "d", "e(bot)"}
	if fourE {
		ekeys = []string{"a(top)", "b", "c", "d(bot)"}
	}
	for i, l := range ekeys {
		out.Plot("t", "sl", l, S[i], -1)
	}

	// lrm
	out.Splot("liquid retention model")
	for i, l := range ekeys {
		out.Plot("pc", "sl", l, S[i], -1)
	}

	// save
	plt.SetForPng(1.2, 600, 150)
	out.Draw("/tmp", "fig_"+fnkey+".png", false, nil)
}
