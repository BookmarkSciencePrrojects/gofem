// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"testing"

	"github.com/cpmech/gofem/ana"
	"github.com/cpmech/gofem/ele/solid"
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/tests"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
)

func Test_sigini01(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("sigini01. zero displacements. initial stresses")

	// fem
	main := fem.NewMain("data/sigini01.sim", "", true, false, false, false, chk.Verbose, 0)

	// set stage
	err := main.SetStage(0)
	if err != nil {
		tst.Errorf("SetStage failed:\n%v", err)
		return
	}

	// initialise solution vectors
	err = main.ZeroStage(0, true)
	if err != nil {
		tst.Errorf("ZeroStage failed:\n%v", err)
		return
	}

	// domain
	dom := main.Domains[0]

	// check displacements
	tolu := 1e-16
	for _, n := range dom.Nodes {
		eqx := n.GetEq("ux")
		eqy := n.GetEq("uy")
		u := []float64{dom.Sol.Y[eqx], dom.Sol.Y[eqy]}
		chk.Vector(tst, "u", tolu, u, nil)
	}

	// analytical solution
	qnV, qnH := -100.0, -50.0
	ν := 0.25
	σx, σy := qnH, qnV
	σz := ν * (σx + σy)
	σref := []float64{σx, σy, σz, 0}

	// check stresses
	e := dom.Elems[0].(*solid.Solid)
	tols := 1e-13
	for idx, _ := range e.IpsElem {
		σ := e.States[idx].Sig
		io.Pforan("σ = %v\n", σ)
		chk.Vector(tst, "σ", tols, σ, σref)
	}
}

func Test_sigini02(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("sigini02. initial stresses. run simulation")

	// fem
	main := fem.NewMain("data/sigini02.sim", "", true, false, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed\n%v", err)
		return
	}

	// domain
	dom := main.Domains[0]

	// external force
	e := dom.Elems[0].(*solid.Solid)
	la.PrintMat("K", e.K, "%10.2f", false)

	// solution
	var sol ana.CteStressPstrain
	sol.Init(fun.Prms{
		&fun.Prm{N: "qnH0", V: -20},
		&fun.Prm{N: "qnV0", V: -20},
		&fun.Prm{N: "qnH", V: -50},
		&fun.Prm{N: "qnV", V: -100},
	})

	// check displacements
	t := dom.Sol.T
	tolu := 1e-16
	for _, n := range dom.Nodes {
		eqx := n.GetEq("ux")
		eqy := n.GetEq("uy")
		u := []float64{dom.Sol.Y[eqx], dom.Sol.Y[eqy]}
		sol.CheckDispl(tst, t, u, n.Vert.C, tolu)
	}

	// check stresses
	tols := 1e-13
	for idx, ip := range e.IpsElem {
		x := e.Cell.Shp.IpRealCoords(e.X, ip)
		σ := e.States[idx].Sig
		io.Pforan("σ = %v\n", σ)
		sol.CheckStress(tst, t, σ, x, tols)
	}
}

func Test_square01(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("square01. ini stress free square")

	// fem
	main := fem.NewMain("data/square01.sim", "", true, false, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed\n%v", err)
		return
	}

	// domain
	dom := main.Domains[0]

	// external force
	e := dom.Elems[0].(*solid.Solid)
	la.PrintMat("K", e.K, "%10.2f", false)

	// solution
	var sol ana.CteStressPstrain
	sol.Init(fun.Prms{
		&fun.Prm{N: "qnH", V: -50},
		&fun.Prm{N: "qnV", V: -100},
	})

	// check displacements
	t := dom.Sol.T
	tolu := 1e-16
	for _, n := range dom.Nodes {
		eqx := n.GetEq("ux")
		eqy := n.GetEq("uy")
		u := []float64{dom.Sol.Y[eqx], dom.Sol.Y[eqy]}
		io.Pfyel("u = %v\n", u)
		sol.CheckDispl(tst, t, u, n.Vert.C, tolu)
	}

	// check stresses
	tols := 1e-13
	for idx, ip := range e.IpsElem {
		x := e.Cell.Shp.IpRealCoords(e.X, ip)
		σ := e.States[idx].Sig
		io.Pforan("σ = %v\n", σ)
		sol.CheckStress(tst, t, σ, x, tols)
	}
}

func Test_selfweight01(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("selfweight01. self-weight")

	// fem
	main := fem.NewMain("data/selfweight01.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed\n%v", err)
		return
	}

	// displacement @ top
	dom := main.Domains[0]
	nod := dom.Vid2node[0]
	eqy := nod.GetEq("uy")
	uy := dom.Sol.Y[eqy]
	uy_cor := -8.737017006803450E-05
	io.Pforan("uy @ top = %v (%v)\n", uy, uy_cor)
	chk.Scalar(tst, "uy @ top", 1e-11, uy, uy_cor)

	// check
	if true {
		skipK := true
		tolK := 1e-17
		tolu := 1e-11
		tols := 1e-5
		tests.CompareResults(tst, "data/selfweight01.sim", "cmp/singleq9grav.cmp", "", tolK, tolu, tols, skipK, chk.Verbose, nil)
	}
}

func Test_selfweight02(tst *testing.T) {

	//tests.Verbose()
	chk.PrintTitle("selfweight02. self-weight")

	// fem
	main := fem.NewMain("data/selfweight02.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := main.Run()
	if err != nil {
		tst.Errorf("Run failed\n%v", err)
		return
	}

	// domain and element
	dom := main.Domains[0]
	ele := dom.Elems[0].(*solid.Solid)

	// solution
	var sol ana.ConfinedSelfWeight
	sol.Init(fun.Prms{
		&fun.Prm{N: "E", V: 1e3},
		&fun.Prm{N: "nu", V: 0.25},
		&fun.Prm{N: "rho", V: 2.0},
		&fun.Prm{N: "g", V: 10.0},
		&fun.Prm{N: "h", V: 1.0},
		&fun.Prm{N: "w", V: 1.0},
	})

	// check displacements
	t := 1.0
	tolu := 1e-16
	for _, n := range dom.Nodes {
		eqx := n.GetEq("ux")
		eqy := n.GetEq("uy")
		eqz := n.GetEq("uz")
		u := []float64{dom.Sol.Y[eqx], dom.Sol.Y[eqy], dom.Sol.Y[eqz]}
		io.Pfyel("x=%v u=%v\n", n.Vert.C, u)
		sol.CheckDispl(tst, t, u, n.Vert.C, tolu)
	}

	// check stresses
	tols := 1e-13
	for idx, ip := range ele.IpsElem {
		x := ele.Cell.Shp.IpRealCoords(ele.X, ip)
		σ := ele.States[idx].Sig
		//io.Pforan("\nσ = %v\n", σ)
		sol.CheckStress(tst, t, σ, x, tols)
	}
}
