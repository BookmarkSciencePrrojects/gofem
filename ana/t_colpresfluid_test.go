// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package ana

import (
	"math"
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func Test_colpresfluid01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("colpresfluid01. pressure on fluid along column")

	R0 := 1.0
	p0 := 0.0
	//C := 4.545454545454545e-07
	C := 1e-2
	H := 10.0
	g := 10.0

	var col ColumnFluidPressure
	col.Init(R0, p0, C, g, H, true)

	tol := 1e-4
	np := 11
	dz := H / float64(np-1)
	io.PfWhite("%8s%14s%14s%14s%14s%23s\n", "z", "p_ana", "R_ana", "p_num", "R_num", "err_p")
	for i := 0; i < np; i++ {
		z := H - float64(i)*dz
		p_ana, R_ana := col.Calc(z)
		p_num, R_num := col.CalcNum(z)
		err_p := math.Abs(p_ana - p_num)
		io.Pf("%8.4f%14.8f%14.8f%14.8f%14.8f%23.15e\n", z, p_ana, R_ana, p_num, R_num, err_p)
		chk.AnaNum(tst, "p", tol, p_ana, p_num, false)
	}

	np = 101
	Z := utl.LinSpace(0, H, np)
	P_ana := make([]float64, np)
	R_ana := make([]float64, np)
	P_num := make([]float64, np)
	R_num := make([]float64, np)
	for i, z := range Z {
		P_ana[i], R_ana[i] = col.Calc(z)
		P_num[i], R_num[i] = col.CalcNum(z)
	}

	pMaxLin := R0 * g * H

	if chk.Verbose {

		plt.Subplot(2, 1, 1)
		plt.Plot(P_num, Z, &plt.A{C: "r", Ls: "-", L: "num"})
		plt.Plot(P_ana, Z, &plt.A{C: "b", Ls: ".", L: "ana", Me: 20})
		plt.Plot([]float64{p0, pMaxLin}, []float64{H, 0}, &plt.A{C: "k", Ls: "--"})
		plt.Gll("$p$", "$z$", nil)

		plt.Subplot(2, 1, 2)
		plt.Plot(R_num, Z, &plt.A{C: "r", Ls: "-", L: "num"})
		plt.Plot(R_ana, Z, &plt.A{C: "b", Ls: ".", L: "ana", Me: 20})
		plt.Plot([]float64{R0, R0 + C*pMaxLin}, []float64{H, 0}, &plt.A{C: "k", Ls: "--"})
		plt.Gll("$\\rho$", "$z$", nil)

		plt.SaveD("/tmp/gofem", "fig_colpresfluid01.eps")
	}
}
