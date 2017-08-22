// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package tests implements structures and functions to test elements and FE simulations
package tests

import (
	"encoding/json"
	"math"
	"testing"

	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/ele/solid"
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

// Iteration holds results from iterations
type Iteration struct {
	It     int     // iteration number
	ResRel float64 // relative residual
	Resid  float64 // absolute residual
}

// Results holds numerical results
type Results struct {
	Status     string        // status message
	LoadFactor float64       // load factor
	Iterations []Iteration   // iterations data
	Kmats      [][][]float64 // [nele][nu][nu] all stiffness matrices
	Disp       [][]float64   // [nnod][ndim] displacements at nodes
	DispMult   float64       // displacements multiplier
	Note       string        // note about number of integration points
	Stresses   [][][]float64 // [nele][nip][nsig] all stresses @ all ips
	// Stresses: examples:
	//   2D solid: [sx, sy, sxy, sz]
	//   3D beam:  [M22, M11, Mtt, Vs, Vr]
}

// ResultsSet is a set of comparison results
type ResultsSet []*Results

// CompareResults performs comparison of results (gofem versus .cmp files)
func CompareResults(tst *testing.T, simfilepath, cmpfname, alias string, tolK, tolu, tols float64, skipK, verbose bool, extraAfterSetStage func(dom *fem.Domain)) {

	// flag
	do_check_stresses := true

	// FEM structure
	main := fem.NewMain(simfilepath, alias, false, false, true, false, verbose, 0)

	// set stage
	err := main.SetStage(0)
	if err != nil {
		chk.Panic("cannot set stage:\n%v", err)
	}

	// extra settings
	if extraAfterSetStage != nil {
		extraAfterSetStage(main.Domains[0])
	}

	// zero solution
	err = main.ZeroStage(0, true)
	if err != nil {
		chk.Panic("cannot zero stage data:\n%v", err)
	}

	// read file with comparison results
	buf, err := io.ReadFile(cmpfname)
	if err != nil {
		tst.Errorf("CompareResults: ReadFile failed:%v\n", err)
		return
	}

	// unmarshal json
	var cmp_set ResultsSet
	err = json.Unmarshal(buf, &cmp_set)
	if err != nil {
		tst.Errorf("CompareResults: Unmarshal failed\n")
		return
	}

	// run comparisons
	dom := main.Domains[0]
	dmult := 1.0
	for idx, cmp := range cmp_set {

		// displacements multiplier
		if idx == 0 && math.Abs(cmp.DispMult) > 1e-10 {
			dmult = cmp.DispMult
		}

		// time index
		tidx := idx + 1
		if verbose {
			io.PfYel("\n\ntidx = %d . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .\n", tidx)
		}

		// load gofem results
		err = dom.Read(main.Summary, tidx, 0, true)
		if err != nil {
			chk.Panic("cannot read 'gofem' results:\n%v", err)
		}
		if verbose {
			io.Pfyel("time = %v\n", dom.Sol.T)
		}

		// check K matrices
		if !skipK {
			if verbose {
				io.Pfgreen(". . . checking K matrices . . .\n")
			}
			for eid, Ksg := range cmp.Kmats {
				if e, ok := dom.Elems[eid].(*solid.Solid); ok {
					err = e.AddToKb(dom.Kb, dom.Sol, true)
					if err != nil {
						chk.Panic("CompareResults: AddToKb failed\n")
					}
					chk.Deep2(tst, io.Sf("K%d", eid), tolK, e.K, Ksg)
				}
				if e, ok := dom.Elems[eid].(*solid.Beam); ok {
					err = e.AddToKb(dom.Kb, dom.Sol, true)
					if err != nil {
						chk.Panic("CompareResults: AddToKb failed\n")
					}
					chk.Deep2(tst, io.Sf("K%d", eid), tolK, e.K, Ksg)
				}
				if e, ok := dom.Elems[eid].(*solid.Rod); ok {
					err = e.AddToKb(dom.Kb, dom.Sol, true)
					if err != nil {
						chk.Panic("CompareResults: AddToKb failed\n")
					}
					chk.Deep2(tst, io.Sf("K%d", eid), tolK, e.K, Ksg)
				}
				if e, ok := dom.Elems[eid].(*solid.ElastRod); ok {
					err = e.AddToKb(dom.Kb, dom.Sol, true)
					if err != nil {
						chk.Panic("CompareResults: AddToKb failed\n")
					}
					chk.Deep2(tst, io.Sf("K%d", eid), tolK, e.K, Ksg)
				}
			}
		}

		// check displacements
		if verbose {
			io.Pfgreen(". . . checking displacements . . .\n")
		}
		for nid, usg := range cmp.Disp {
			ix := dom.Vid2node[nid].Dofs[0].Eq
			iy := dom.Vid2node[nid].Dofs[1].Eq
			chk.AnaNum(tst, "ux", tolu, dom.Sol.Y[ix], usg[0]*dmult, verbose)
			chk.AnaNum(tst, "uy", tolu, dom.Sol.Y[iy], usg[1]*dmult, verbose)
			if len(usg) > 2 {
				iz := dom.Vid2node[nid].Dofs[2].Eq
				chk.AnaNum(tst, "uz", tolu, dom.Sol.Y[iz], usg[2]*dmult, verbose)
				for j := 3; j < len(usg); j++ {
					idx := dom.Vid2node[nid].Dofs[j].Eq
					chk.AnaNum(tst, io.Sf("u%d", j), tolu, dom.Sol.Y[idx], usg[j]*dmult, verbose)
				}
			}
		}

		// check stresses
		if do_check_stresses {
			if verbose {
				io.Pfgreen(". . . checking stresses (or moments/shear forces). . .\n")
			}
			for eid, sig := range cmp.Stresses {
				if verbose {
					io.Pforan("eid = %d\n", eid)
				}
				if e, ok := dom.Cid2elem[eid].(*solid.Solid); ok {
					for ip, val := range sig {
						if verbose {
							io.Pfgrey2("ip = %d @ %v\n", ip, e.Cell.Shp.IpRealCoords(e.X, e.IpsElem[ip]))
						}
						σ := e.States[ip].Sig
						if len(val) == 6 {
							chk.AnaNum(tst, "sx ", tols, σ[0], val[0], verbose)
							chk.AnaNum(tst, "sy ", tols, σ[1], val[1], verbose)
						} else {
							chk.AnaNum(tst, "sx ", tols, σ[0], val[0], verbose)
							chk.AnaNum(tst, "sy ", tols, σ[1], val[1], verbose)
							chk.AnaNum(tst, "sxy", tols, σ[3]/math.Sqrt2, val[2], verbose)
							if len(val) > 3 { // sx, sy, sxy, sz
								chk.AnaNum(tst, "sz ", tols, σ[2], val[3], verbose)
							}
						}
					}
				}
				if e, ok := dom.Cid2elem[eid].(*solid.Rod); ok {
					for ip, val := range sig {
						if verbose {
							io.Pfgrey2("ip = %d\n", ip)
						}
						σ := e.States[ip].Sig
						chk.AnaNum(tst, "sig", tols, σ, val[0], verbose)
					}
				}
				if e, ok := dom.Cid2elem[eid].(*solid.ElastRod); ok {
					res := ele.NewIpsMap()
					e.OutIpVals(res, dom.Sol)
					for ip, val := range sig {
						if verbose {
							io.Pfgrey2("ip = %d\n", ip)
						}
						σ := res.Get("sig", 0)
						chk.AnaNum(tst, "sig", tols, σ, val[0], verbose)
					}
				}
				if e, ok := dom.Cid2elem[eid].(*solid.Beam); ok {
					res := ele.NewIpsMap()
					e.OutIpVals(res, dom.Sol)
					for st, val := range sig {
						if verbose {
							io.Pfgrey2("station = %d\n", st)
						}
						if e.Ndim == 2 {
							M22 := res.Get("M22", st)
							chk.AnaNum(tst, "M22", tols, M22, val[0], verbose)
						} else {
							M22, M11, Mtt := res.Get("M22", st), res.Get("M11", st), res.Get("T00", st)
							chk.AnaNum(tst, "M22", tols, M22, val[0], verbose)
							chk.AnaNum(tst, "M11", tols, M11, val[1], verbose)
							chk.AnaNum(tst, "Mtt", tols, Mtt, val[2], verbose)
						}
					}
				}
			}
		}
	}
}
