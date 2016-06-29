// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be fndim in the LICENSE file.

package seepage

func LiqFlowKeys(ndim int) []string {
	// nwl == nl・wl == filter velocity
	if ndim == 2 {
		return []string{"nwlx", "nwly"}
	}
	return []string{"nwlx", "nwly", "nwlz"}
}

func GasFlowKeys(ndim int) []string {
	// nwg == ng・wg == figter vegocity
	if ndim == 2 {
		return []string{"nwgx", "nwgy"}
	}
	return []string{"nwgx", "nwgy", "nwgz"}
}
