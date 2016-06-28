// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

func main() {

	// filename
	filename, fnkey := io.ArgToFilename(0, "a-coarse-elast-d2-q9.sim", ".sim", true)

	// start analysis process
	out.Start(filename, 0, 0)

	// define entities
	out.Define("A B C D E", out.N{-1, -2, -3, -4, -5})                       // top => bottom
	out.Define("a b c d e", out.P{{18, 8}, {8, 8}, {4, 8}, {30, 8}, {0, 0}}) // top => bottom

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
	out.Splot("t-pl", "liquid pressure")
	for i, l := range []string{"A", "B", "C", "D", "E"} {
		out.Plot("t", "pl", l, S[i], -1)
	}

	// uy
	out.Splot("t-uy", "displacements")
	for i, l := range []string{"A", "B", "C", "D", "E"} {
		out.Plot("t", "uy", l, S[i], -1)
	}

	out.Splot("t-sl", "liquid saturation")
	for i, l := range []string{"a", "b", "c", "d", "e"} {
		out.Plot("t", "sl", l, S[i], -1)
	}

	out.Splot("t-sy", "stresses")
	for i, l := range []string{"a", "b", "c", "d", "e"} {
		out.Plot("t", "sy", l, S[i], -1)
	}

	// show
	plt.SetForPng(1, 500, 200)
	out.Draw("/tmp", "up_indentation2d_unsat_"+fnkey+".png", -1, -1, false, nil)
}
