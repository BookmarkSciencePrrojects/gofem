// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"math"

	"github.com/cpmech/gofem/ana"
	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/ele/solid"
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func main() {

	// catch errors
	defer func() {
		if err := recover(); err != nil {
			io.PfRed("ERROR: %v\n", err)
		}
	}()

	// filename
	filename, fnkey := io.ArgToFilename(0, "spo751", ".sim", true)

	// constants
	nidx := 20 // selected node at outer surface
	didx := 0  // selected  dof index for plot
	nels := 4  // number of elements
	nips := 4  // number of ips

	// selected P values for stress plot
	Psel := []float64{100, 140, 180, 190}
	tolPsel := 2.0    // tolerance to compare P
	GPa2MPa := 1000.0 // conversion factor

	// input data
	Pcen := 200.0         // [Mpa]
	a, b := 100.0, 200.0  // [mm], [mm]
	E, ν := 210000.0, 0.3 // [MPa], [-]
	σy := 240.0           // [MPa]

	// analytical solution
	var sol ana.PressCylin
	sol.Init([]*fun.P{
		&fun.P{N: "a", V: a}, &fun.P{N: "b", V: b},
		&fun.P{N: "E", V: E}, &fun.P{N: "ν", V: ν},
		&fun.P{N: "σy", V: σy},
	})
	np := 41
	P_ana, Ub_ana := sol.CalcPressDisp(np)
	R_ana, Sr_ana, St_ana := sol.CalcStresses(Psel, np)

	// fem
	analysis := fem.NewMain(filename, "", false, false, true, false, true, 0)
	err := analysis.SetStage(0)
	if err != nil {
		chk.Panic("SetStage failed:\n%v", err)
	}
	err = analysis.ZeroStage(0, true)
	if err != nil {
		chk.Panic("ZeroStage failed:\n%v", err)
	}
	dom := analysis.Domains[0]
	sum := analysis.Summary

	// gofem results
	nto := len(sum.OutTimes)
	P := make([]float64, nto)
	Ub := make([]float64, nto)
	R := utl.Deep3alloc(len(Psel), nels, nips)
	Sr := utl.Deep3alloc(len(Psel), nels, nips)
	St := utl.Deep3alloc(len(Psel), nels, nips)
	i := 0
	for tidx, t := range sum.OutTimes {

		// read results from file
		err = dom.Read(sum, tidx, 0, true)
		if err != nil {
			chk.Panic("cannot read solution\n%v", err)
		}

		// collect results for load versus displacement plot
		nod := dom.Nodes[nidx]
		eq := nod.Dofs[didx].Eq
		P[tidx] = t * Pcen
		Ub[tidx] = dom.Sol.Y[eq]

		// stresses
		if isPsel(Psel, P[tidx], tolPsel) {
			for j, ee := range dom.ElemIntvars {
				e := ee.(*solid.Solid)
				X := e.OutIpCoords()
				M := ele.NewIpsMap()
				e.OutIpVals(M, dom.Sol)
				for k := 0; k < nips; k++ {
					x, y := X[k][0], X[k][1]
					sx := M.Get("sx", k) * GPa2MPa
					sy := M.Get("sy", k) * GPa2MPa
					sxy := M.Get("sxy", k) * GPa2MPa / math.Sqrt2
					R[i][j][k], Sr[i][j][k], St[i][j][k], _ = ana.PolarStresses(x, y, sx, sy, sxy)
				}
			}
			i++
		}
	}

	// auxiliary data for plotting stresses
	colors := []string{"r", "m", "g", "k", "y", "c", "r", "m"}
	markers1 := []string{"o", "s", "x", ".", "^", "*", "o", "s"}
	markers2 := []string{"+", "+", "+", "+", "+", "+", "+", "+"}

	// plot load displacements
	if true {
		//if false {
		plt.Plot(Ub_ana, P_ana, &plt.A{C: "b", Ms: 2, L: "solution"})
		plt.Plot(Ub, P, &plt.A{C: "r", M: ".", Ls: "--", L: "fem: outer"})
		plt.Gll("$u_x\\;\\mathrm{[mm]}$", "$P\\;\\mathrm{[MPa]}$", nil)
		plt.Save("/tmp", io.Sf("gofem_%s_disp", fnkey))
	}

	// plot radial stresses
	if true {
		//if false {
		plt.Reset(false, nil)
		for i, Pval := range Psel {
			plt.Plot(R_ana, Sr_ana[i], &plt.A{C: "b"})
			for k := 0; k < nips; k++ {
				for j := 0; j < nels; j++ {
					args := &plt.A{C: colors[i], M: markers1[i]}
					if k > 1 {
						args.C = "k"
						args.M = markers2[i]
					}
					if k == 0 && j == 0 {
						args.L = io.Sf("P=%g", Pval)
					}
					plt.PlotOne(R[i][j][k], Sr[i][j][k], args)
				}
			}
		}
		plt.Gll("$r\\;\\mathrm{[mm]}$", "$\\sigma_r\\;\\mathrm{[MPa]}$", &plt.A{LegLoc: "lower right"})
		plt.AxisXrange(a, b)
		plt.Save("/tmp", io.Sf("gofem_%s_sr", fnkey))
	}

	// plot tangential stresses
	if true {
		//if false {
		plt.Reset(false, nil)
		for i, Pval := range Psel {
			plt.Plot(R_ana, St_ana[i], &plt.A{C: "b"})
			for k := 0; k < nips; k++ {
				for j := 0; j < nels; j++ {
					args := &plt.A{C: colors[i], M: markers1[i]}
					if k > 1 {
						args.C = "k"
						args.M = markers2[i]
					}
					if k == 0 && j == 0 {
						args.L = io.Sf("P=%g", Pval)
					}
					plt.PlotOne(R[i][j][k], St[i][j][k], args)
				}
			}
		}
		plt.Gll("$r\\;\\mathrm{[mm]}$", "$\\sigma_t\\;\\mathrm{[MPa]}$", &plt.A{LegLoc: "upper left"})
		plt.Save("/tmp", io.Sf("gofem_%s_st", fnkey))
	}
}

func isPsel(Psel []float64, p, tol float64) bool {
	for _, pp := range Psel {
		if math.Abs(p-pp) < tol {
			return true
		}
	}
	return false
}
