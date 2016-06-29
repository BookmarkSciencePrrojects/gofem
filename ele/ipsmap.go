// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package ele

// IpsMap defines a map to hold results @ integration points
type IpsMap map[string][]float64

// NewIpsMap returns a new IpsMap
func NewIpsMap() *IpsMap {
	var M IpsMap
	M = make(map[string][]float64)
	return &M
}

// Set sets item in map by key and ip-index. The slice is resized with nip in case it's empty
//  Input:
//   idx -- index of integration point
//   nip -- number of integration points (to resize if necessary)
//   val -- value of 'key' @ integration point 'idx'
func (o *IpsMap) Set(key string, idx, nip int, val float64) {
	if slice, ok := (*o)[key]; ok {
		slice[idx] = val
		return
	}
	slice := make([]float64, nip)
	slice[idx] = val
	(*o)[key] = slice
}

// Get returns item corresponding to 'key' and integration point 'idx'
//  Note: this function returns 0 if 'key' is not found. It also does not check for out-of-bound errors
func (o *IpsMap) Get(key string, idx int) float64 {
	if slice, ok := (*o)[key]; ok {
		return slice[idx]
	}
	return 0
}
