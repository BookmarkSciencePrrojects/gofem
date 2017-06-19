// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package ele

import "github.com/cpmech/gosl/fun/dbf"

// NaturalBc holds information on natural boundary conditioins such as
// distributed loads or fluxes acting on surfaces
type NaturalBc struct {
	Key     string // key such as qn, qn0, ql, seepH, seepP, etc...
	IdxFace int    // local index of face
	Fcn     dbf.T  // function callback
	Extra   string // extra information
}
