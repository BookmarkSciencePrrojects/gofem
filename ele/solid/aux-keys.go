// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be fndim in the LICENSE file.

package solid

func StressKeys(ndim int) []string {
	if ndim == 2 {
		return []string{"sx", "sy", "sz", "sxy"}
	}
	return []string{"sx", "sy", "sz", "sxy", "syz", "szx"}
}

// Ivs2sigmas converts ivs map to σ values [nsig]
//  σ -- [ndim] stresses
//  i -- index of integration point
func Ivs2sigmas(σ []float64, i int, ivs map[string][]float64) {
	for key, vals := range ivs {
		switch key {
		case "sx":
			σ[0] = vals[i]
		case "sy":
			σ[1] = vals[i]
		case "sz":
			σ[2] = vals[i]
		case "sxy":
			σ[3] = vals[i]
		case "syz":
			if len(σ) > 4 {
				σ[4] = vals[i]
			}
		case "szx":
			if len(σ) > 5 {
				σ[5] = vals[i]
			}
		}
	}
}
