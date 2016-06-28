// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package por

// GetPathCycle sets a path with cycles of drying-wetting
//  pc0   -- initial capillary pressure
//  P     -- example: {10, 5, 20, 0}
//  np    -- number of points in each branch
func GetPathCycle(pc0 float64, P []float64, np int) (Pc []float64) {
	if np < 2 {
		np = 2
	}
	Pc = []float64{pc0}
	k := 1
	for _, pc := range P {
		δpc := (pc - Pc[k-1]) / float64(np-1)
		for i := 1; i < np; i++ {
			Pc = append(Pc, Pc[k-1]+δpc)
			k += 1
		}
	}
	return Pc
}
