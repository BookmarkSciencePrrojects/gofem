// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"strings"

	"github.com/cpmech/gofem/mdl/lrm"
	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func main() {

	// filename
	filename, fnkey := io.ArgToFilename(0, "upp_3mcol_desic16e.sim", ".sim", true)
	has16e := strings.HasSuffix(fnkey, "desic16e")
	upForm := strings.HasPrefix(fnkey, "up_3mcol")

	// start analysis process
	out.Start(filename, 0, 0)

	// labels
	labelsN := []string{"A(top)", "B", "C(mid)", "D", "E(bot)"}
	labelsE := []string{"a(top)", "b", "c", "d(bot)"}
	if has16e {
		labelsE = []string{"a(top)", "b", "c(mid)", "d", "e(bot)"}
	}

	// define entities
	out.Define("A(top) B C(mid) D E(bot)", out.N{-5, -4, -3, -2, -1})
	if has16e {
		out.Define("a(top) b c(mid) d e(bot)", out.P{{15, 8}, {12, 4}, {8, 4}, {4, 4}, {0, 0}})
	} else {
		out.Define("a(top) b c d(bot)", out.P{{3, 8}, {2, 4}, {1, 4}, {0, 0}})
	}
	out.Define("left-side", out.Along{{0, 0}, {0, 3}})

	// load results
	out.LoadResults(nil)

	// styles
	me := 10
	S := []plt.Fmt{
		plt.Fmt{C: "g", M: "o", Me: me},
		plt.Fmt{C: "b", M: "*", Me: me},
		plt.Fmt{C: "m", M: "x", Me: me},
		plt.Fmt{C: "orange", M: "+", Me: me},
		plt.Fmt{C: "r", M: "^", Me: me},
	}

	// time outputs
	//io.Pforan("out.Times = %v\n", out.Times)
	I, _ := utl.GetITout(out.Times, []float64{0, 340, 700, 1000, 1400}, 0.1)

	// plot LRM only
	onlyLRM := false

	// all plots
	if !onlyLRM {

		// pl versus y
		out.Splot("pl-y", "liquid pressure along column")
		for _, i := range I {
			t := out.Times[i]
			out.Plot("pl", "y", "left-side", plt.Fmt{L: io.Sf("t=%g", t)}, i)
		}

		// pg versus y
		if upForm {
			out.Splot("empty1", "")
		} else {
			out.Splot("pg-y", "gas pressure along column")
			for _, i := range I {
				t := out.Times[i]
				out.Plot("pg", "y", "left-side", plt.Fmt{L: io.Sf("t=%g", t)}, i)
			}
		}

		// pl versus t
		out.Splot("pl-t", "liquid pressure in time")
		for i, l := range labelsN {
			out.Plot("t", "pl", l, S[i], -1)
		}

		// pg versus t
		if upForm {
			out.Splot("empty2", "")
		} else {
			out.Splot("t-pg", "gas pressure in time")
			for i, l := range labelsN {
				out.Plot("t", "pg", l, S[i], -1)
			}
		}

		// uy versus t
		out.Splot("t-uy", "displacements in time")
		for i, l := range labelsN {
			out.Plot("t", "uy", l, S[i], -1)
		}

		// sl versus t
		out.Splot("t-sl", "liquid saturation in time")
		for i, l := range labelsE {
			out.Plot("t", "sl", l, S[i], -1)
		}

		// nwly versus time
		out.Splot("t-nwly", "liquid filter velocity versus time")
		for i, l := range labelsE {
			out.Plot("t", "nwly", l, S[i], -1)
		}

		// nwgy versus time
		if upForm {
			out.Splot("empty3", "")
		} else {
			out.Splot("t-nwgy", "gas filter velocity versus time")
			for i, l := range labelsE {
				out.Plot("t", "nwgy", l, S[i], -1)
			}
		}

		// sy versus t
		out.Splot("t-sy", "vertical stress in time")
		for i, l := range labelsE {
			out.Plot("t", "sy", l, S[i], -1)
		}

	}

	// pc versus sl
	out.Splot("lrm", "liquid retention behaviour")
	for i, l := range labelsE {
		out.Plot("pc", "sl", l, S[i], -1)
	}

	// save
	sim := out.Dom.Sim
	if onlyLRM {
		plt.SetForPng(0.75, 400, 200)
	} else {
		plt.SetForPng(1.8, 600, 200)
	}
	out.Draw("/tmp/gofem", "fig_"+fnkey+".png", -1, -1, false, func(id string) {
		if id == "lrm" {
			mat := sim.MatModels.Get("lreten1")
			Lrm := mat.Lrm
			lrm.Plot(Lrm, 0, Lrm.SlMax(), 30, 101, "'k-^', markerfacecolor='white', ms=5, markevery=10", "", "model")
		}
	})
}
