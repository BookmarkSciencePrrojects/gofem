// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package ele

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/la"
)

// BuildCoordsMatrix returns the coordinate matrix of a particular Cell
func BuildCoordsMatrix(cell *inp.Cell, msh *inp.Mesh) (x [][]float64) {
	x = la.MatAlloc(msh.Ndim, len(cell.Verts))
	for i := 0; i < msh.Ndim; i++ {
		for j, v := range cell.Verts {
			x[i][j] = msh.Verts[v].C[i]
		}
	}
	return
}
