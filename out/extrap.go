// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/ele/porous"
	"github.com/cpmech/gofem/ele/seepage"
	"github.com/cpmech/gofem/ele/solid"
	"github.com/cpmech/gofem/shp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/la"
)

func ComputeExtrapolatedValues(extrapKeys []string) {

	// auxiliary
	verts := Dom.Msh.Verts
	cells := Dom.Msh.Cells

	// allocate structures for extrapolation
	nverts := len(verts)
	ExVals = make([]map[string]float64, nverts)
	counts := make([]map[string]float64, nverts)
	for i := 0; i < nverts; i++ {
		ExVals[i] = make(map[string]float64)
		counts[i] = make(map[string]float64)
	}

	// loop over elements
	for _, element := range Dom.Elems {
		if e, ok := element.(ele.CanOutputIps); ok {

			// get shape and integration points from known elements
			var sha *shp.Shape
			var ips []shp.Ipoint
			switch e := element.(type) {
			case *seepage.Liquid:
				sha = e.Cell.Shp
				ips = e.IpsElem
			case *solid.Solid:
				sha = e.Cell.Shp
				ips = e.IpsElem
			case *porous.SolidLiquid:
				sha = e.U.Cell.Shp
				ips = e.U.IpsElem
			case *porous.SolidLiquidGas:
				sha = e.U.Cell.Shp
				ips = e.U.IpsElem
			}
			if sha == nil {
				return // cannot extrapolate; e.g. rjoint, beams
				//chk.Panic("cannot get shape structure from element")
			}

			// compute Extrapolator matrix
			Emat := la.MatAlloc(sha.Nverts, len(ips))
			err := sha.Extrapolator(Emat, ips)
			if err != nil {
				chk.Panic("cannot compute extrapolator matrix: %v", err)
			}

			// get ips data
			allvals := ele.NewIpsMap()
			e.OutIpVals(allvals, Dom.Sol)

			// perform extrapolation
			cell := cells[element.Id()]
			for _, key := range extrapKeys {
				if vals, ok := (*allvals)[key]; ok {
					for i := 0; i < sha.Nverts; i++ {
						v := cell.Verts[i]
						for j := 0; j < len(ips); j++ {
							ExVals[v][key] += Emat[i][j] * vals[j]
						}
					}
				} else {
					chk.Panic("ip does not have key = %s", key)
				}
			}

			// increment counter
			for i := 0; i < sha.Nverts; i++ {
				v := cell.Verts[i]
				for _, key := range extrapKeys {
					counts[v][key] += 1
				}
			}
		}
	}

	// compute average
	for i := 0; i < nverts; i++ {
		for key, cnt := range counts[i] {
			ExVals[i][key] /= cnt
		}
	}
}
