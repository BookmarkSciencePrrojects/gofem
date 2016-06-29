// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package out implements FE simulation output handling for analyses and plotting
package out

import (
	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/ele/solid"
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/gm"
	"github.com/cpmech/gosl/utl"
)

// constants
var (
	TolC = 1e-8 // tolerance to compare x-y-z coordinates
	TolT = 1e-3 // tolerance to compare times
	Ndiv = 20   // bins n-division
)

// ResultsMap maps aliases to points
type ResultsMap map[string]Points

// IpData_t is an auxiliary structure holding element id and integration point coordinates
type IpData_t struct {
	Cid  int                // id of cell/element holding this integration point
	X    []float64          // coordinates of integration point
	Vals map[string]float64 // current (@ time t) values
}

// Global variables
var (

	// data set by Start
	Analysis   *fem.Main          // the fem structure
	Sum        *fem.Summary       // [from Analysis] summary
	Dom        *fem.Domain        // [from Analysis] FE domain
	Ipoints    []*IpData_t        // all integration points. ipid == index in Ipoints
	Cid2ips    [][]int            // [ncells][nip] maps cell id to index in Ipoints
	Ipkey2ips  map[string][]int   // maps ip keys to indices in Ipoints
	Ipkeys     map[string]bool    // all ip keys
	NodBins    gm.Bins            // bins for nodes
	IpsBins    gm.Bins            // bins for integration points
	IpsMin     []float64          // [ndim] {x,y,z}_min among all ips
	IpsMax     []float64          // [ndim] {x,y,z}_max among all ips
	Beams      []*solid.Beam      // beams, if any
	ElemOutIps []ele.CanOutputIps // subset of element that can output IP values

	// defined entities and results loaded by LoadResults
	Planes   map[string]*PlaneData // for points defined on planes. maps aliases to data
	Results  ResultsMap            // maps labels => points
	TimeInds []int                 // selected output indices
	Times    []float64             // selected output times

	// extrapolated values
	Extrap []string             // keys to be extrapolated; e.g. []string{"nwlx", "nwly"}
	ExVals []map[string]float64 // [nverts][nkeys] extrapolated values

	// subplots
	Splots []*SplotDat // all subplots
	Csplot *SplotDat   // current subplot
)

// Start starts handling of results given a simulation input file
func Start(simfnpath string, stageIdx, regionIdx int) {

	// fem structure
	Analysis = fem.NewMain(simfnpath, "", false, false, true, false, false, 0)
	Dom = Analysis.Domains[regionIdx]
	Sum = Analysis.Summary

	// set stage
	err := Analysis.SetStage(stageIdx)
	if err != nil {
		chk.Panic("cannot set stage:\n%v", err)
	}

	// initialise solution vectors
	err = Analysis.ZeroStage(stageIdx, true)
	if err != nil {
		chk.Panic("cannot initialise solution vectors:\n%v", err)
	}

	// clear previous data
	Ipoints = make([]*IpData_t, 0)
	Cid2ips = make([][]int, len(Dom.Msh.Cells))
	Ipkey2ips = make(map[string][]int)
	Ipkeys = make(map[string]bool)
	Planes = make(map[string]*PlaneData)
	Results = make(map[string]Points)
	TimeInds = make([]int, 0)
	Times = make([]float64, 0)
	Splots = make([]*SplotDat, 0)
	Beams = make([]*solid.Beam, 0)
	ElemOutIps = make([]ele.CanOutputIps, 0)

	// bins
	m := Dom.Msh
	δ := TolC * 2
	xi := []float64{m.Xmin - δ, m.Ymin - δ}
	xf := []float64{m.Xmax + δ, m.Ymax + δ}
	if m.Ndim == 3 {
		xi = append(xi, m.Zmin-δ)
		xf = append(xf, m.Zmax+δ)
	}
	err = NodBins.Init(xi, xf, Ndiv)
	if err != nil {
		chk.Panic("cannot initialise bins for nodes: %v", err)
	}
	err = IpsBins.Init(xi, xf, Ndiv)
	if err != nil {
		chk.Panic("cannot initialise bins for integration points: %v", err)
	}

	// add nodes to bins
	for _, nod := range Dom.Nodes {
		err := NodBins.Append(nod.Vert.C, nod.Vert.Id)
		if err != nil {
			return
		}
	}

	// to find limits
	ndim := Dom.Msh.Ndim
	IpsMin = make([]float64, ndim)
	IpsMax = make([]float64, ndim)
	first := true

	// for all cells/elements
	for cid, element := range Dom.Cid2elem {
		if element == nil {
			continue
		}

		// add integration points to slice of ips and to bins
		if e, ok := element.(ele.CanOutputIps); ok {
			coords := e.OutIpCoords()
			keys := e.OutIpKeys()
			nip := len(coords)
			ids := make([]int, nip)
			for i := 0; i < nip; i++ {

				// add to bins
				ipid := len(Ipoints)
				ids[i] = ipid
				Ipoints = append(Ipoints, &IpData_t{cid, coords[i], make(map[string]float64)})
				err = IpsBins.Append(coords[i], ipid)
				if err != nil {
					chk.Panic("cannot append to bins of integration points: %v", err)
				}

				// set auxiliary map
				for _, key := range keys {
					utl.StrIntsMapAppend(&Ipkey2ips, key, ipid)
					Ipkeys[key] = true
				}

				// limits
				if first {
					for j := 0; j < ndim; j++ {
						IpsMin[j] = coords[i][j]
						IpsMax[j] = coords[i][j]
					}
					first = false
				} else {
					for j := 0; j < ndim; j++ {
						IpsMin[j] = utl.Min(IpsMin[j], coords[i][j])
						IpsMax[j] = utl.Max(IpsMax[j], coords[i][j])
					}
				}
			}
			Cid2ips[cid] = ids
			ElemOutIps = append(ElemOutIps, e)
		}

		// find beams
		if beam, ok := element.(*solid.Beam); ok {
			Beams = append(Beams, beam)
		}
	}
}
