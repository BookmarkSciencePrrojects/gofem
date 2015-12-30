// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"github.com/cpmech/gofem/fem"
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
	for _, ele := range Dom.Elems {
		if e, ok := ele.(fem.ElemOutIps); ok {

			// get shape and integration points from known elements
			var sha *shp.Shape
			var ips []shp.Ipoint
			switch e := ele.(type) {
			case *fem.ElemP:
				sha = e.Cell.Shp
				ips = e.IpsElem
			case *fem.ElemU:
				sha = e.Cell.Shp
				ips = e.IpsElem
			case *fem.ElemUP:
				sha = e.U.Cell.Shp
				ips = e.U.IpsElem
			case *fem.ElemUPP:
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
			allvals := fem.NewIpsMap()
			e.OutIpVals(allvals, Dom.Sol)

			// perform extrapolation
			cell := cells[ele.Id()]
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
