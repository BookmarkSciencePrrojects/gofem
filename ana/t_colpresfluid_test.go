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
	col.Init(R0, p0, C, g, H)

	tol := 1e-8
	np := 11
	dz := H / float64(np-1)
	io.PfWhite("%8s%14s%14s%14s%14s%23s\n", "z", "pAna", "Rana", "pNum", "Rnum", "errp")
	for i := 0; i < np; i++ {
		z := H - float64(i)*dz
		pAna, Rana := col.Calc(z)
		pNum, Rnum := col.CalcNum(z)
		errp := math.Abs(pAna - pNum)
		io.Pf("%8.4f%14.8f%14.8f%14.8f%14.8f%23.15e\n", z, pAna, Rana, pNum, Rnum, errp)
		chk.AnaNum(tst, "p", tol, pAna, pNum, false)
	}

	np = 101
	Z := utl.LinSpace(0, H, np)
	Pana := make([]float64, np)
	Rana := make([]float64, np)
	Pnum := make([]float64, np)
	Rnum := make([]float64, np)
	for i, z := range Z {
		Pana[i], Rana[i] = col.Calc(z)
		Pnum[i], Rnum[i] = col.CalcNum(z)
	}

	pMaxLin := R0 * g * H

	if chk.Verbose {

		plt.Reset(false, nil)
		plt.Subplot(2, 1, 1)
		plt.Plot(Pnum, Z, &plt.A{C: "r", Ls: "-", L: "num"})
		plt.Plot(Pana, Z, &plt.A{C: "b", Ls: ".", L: "ana", Me: 20})
		plt.Plot([]float64{p0, pMaxLin}, []float64{H, 0}, &plt.A{C: "k", Ls: "--"})
		plt.Gll("$p$", "$z$", nil)

		plt.Subplot(2, 1, 2)
		plt.Plot(Rnum, Z, &plt.A{C: "r", Ls: "-", L: "num"})
		plt.Plot(Rana, Z, &plt.A{C: "b", Ls: ".", L: "ana", Me: 20})
		plt.Plot([]float64{R0, R0 + C*pMaxLin}, []float64{H, 0}, &plt.A{C: "k", Ls: "--"})
		plt.Gll("$\\rho$", "$z$", nil)

		plt.Save("/tmp/gofem", "fig_colpresfluid01")
	}
}
