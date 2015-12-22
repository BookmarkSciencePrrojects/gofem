// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"github.com/cpmech/gofem/mdl/lrm"
	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func main() {

	// filename
	filename, fnkey := io.ArgToFilename(0, "upp_3mcol_desic.sim", ".sim", true)

	// start analysis process
	out.Start(filename, 0, 0)

	// define entities
	out.Define("A(top) B C(mid) D E(bot)", out.N{-5, -4, -3, -2, -1})
	out.Define("a(top) b c d(bot)", out.P{{3, 1}, {2, 8}, {1, 8}, {0, 8}})
	out.Define("left-side", out.Along{{0, 0}, {0, 3}})

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

	// time outputs
	I, _ := utl.GetITout(out.Times, []float64{0, 500, 1000, 1500, 2000}, 0.1)

	if false {

		// pl versus y
		out.Splot("liquid pressure along column")
		for _, i := range I {
			t := out.Times[i]
			out.Plot("pl", "y", "left-side", plt.Fmt{L: io.Sf("t=%g", t)}, i)
		}

		// pg versus y
		out.Splot("gas pressure along column")
		for _, i := range I {
			t := out.Times[i]
			out.Plot("pg", "y", "left-side", plt.Fmt{L: io.Sf("t=%g", t)}, i)
		}

		// pl versus t
		out.Splot("liquid pressure in time")
		for i, l := range []string{"A(top)", "B", "C(mid)", "D", "E(bot)"} {
			out.Plot("t", "pl", l, S[i], -1)
		}

		// pg versus t
		out.Splot("gas pressure in time")
		for i, l := range []string{"A(top)", "B", "C(mid)", "D", "E(bot)"} {
			out.Plot("t", "pg", l, S[i], -1)
		}

		// uy versus t
		out.Splot("displacements in time")
		for i, l := range []string{"A(top)", "B", "C(mid)", "D", "E(bot)"} {
			out.Plot("t", "uy", l, S[i], -1)
		}

		// sl versus t
		out.Splot("liquid saturation in time")
		for i, l := range []string{"a(top)", "b", "c", "d(bot)"} {
			out.Plot("t", "sl", l, S[i], -1)
		}

		// nwly versus time
		out.Splot("liquid filter velocity versus time")
		for i, l := range []string{"a(top)", "b", "c", "d(bot)"} {
			out.Plot("t", "nwly", l, S[i], -1)
		}

		// nwgy versus time
		out.Splot("gas filter velocity versus time")
		for i, l := range []string{"a(top)", "b", "c", "d(bot)"} {
			out.Plot("t", "nwgy", l, S[i], -1)
		}

		// sy versus t
		out.Splot("vertical stress in time")
		for i, l := range []string{"a(top)", "b", "c", "d(bot)"} {
			out.Plot("t", "sy", l, S[i], -1)
		}

	}

	// pc versus sl
	out.Splot("liquid retention behaviour")
	for i, l := range []string{"a(top)", "b", "c", "d(bot)"} {
		out.Plot("pc", "sl", l, S[i], -1)
	}
	//out.Plot("pc", "sl", "d(bot)", S[0], -1)

	// save
	sim := out.Dom.Sim
	//plt.SetForPng(1.8, 800, 150)
	plt.SetForPng(0.75, 400, 200)
	out.Draw("/tmp/gofem", "fig_"+fnkey+".png", false, func(i, j, n int) {
		//if i == 5 && j == 2 {
		mat := sim.MatModels.Get("lreten1")
		Lrm := mat.Lrm
		lrm.Plot(Lrm, 0, Lrm.SlMax(), 30, 101, "'k-^', markerfacecolor='white', ms=5, markevery=10", "", "model")
		//}
	})
}
