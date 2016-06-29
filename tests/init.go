// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package tests

import (
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func init() {
	io.Verbose = false
}

func Verbose() {
	io.Verbose = true
	chk.Verbose = true
}

func GetNidsEqs(dom *fem.Domain) (nids, eqs []int) {
	for _, nod := range dom.Nodes {
		nids = append(nids, nod.Vert.Id)
		for _, dof := range nod.Dofs {
			eqs = append(eqs, dof.Eq)
		}
	}
	return
}
