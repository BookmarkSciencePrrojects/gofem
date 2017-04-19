// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"encoding/json"

	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

func main() {

	// filename
	filename, fnkey := io.ArgToFilename(0, "sg114", ".sim", true)

	// results
	out.Start(filename, 0, 0)
	out.Define("tip", out.N{17})
	out.LoadResults(nil)

	// plot FEM results
	out.Plot("t", "uy", "tip", &plt.A{C: "k", Ls: "-", M: "*", L: "gofem"}, -1)

	// old results
	b, err := io.ReadFile("cmp/sg114gofemold.json")
	if err != nil {
		io.PfRed("cannot read comparison file\n")
		return
	}
	var gofemold struct {
		Time, Uy17 []float64
	}
	err = json.Unmarshal(b, &gofemold)
	if err != nil {
		io.PfRed("cannot unmarshal comparison file\n")
		return
	}

	// mechsys results
	_, res, err := io.ReadTable("cmp/sg114mechsysN17.cmp")
	if err != nil {
		io.PfRed("cannot read mechsys comparison file\n")
		return
	}

	// save
	out.Draw("/tmp", fnkey+".png", -1, -1, false, func(id string) {
		plt.Plot(gofemold.Time, gofemold.Uy17, &plt.A{C: "r", M: ".", L: "gofemOld"})
		plt.Plot(res["Time"], res["uy"], &plt.A{C: "b", M: "+", L: "mechsys"})
	})
}
