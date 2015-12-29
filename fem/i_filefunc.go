// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
)

// IniSetFileFunc sets initial state with values given in a file or by a function
func (o *Domain) IniSetFileFunc(stg *inp.Stage) (err error) {

	// check
	if len(stg.IniFcn.Fcns) != len(stg.IniFcn.Dofs) {
		return chk.Err("number of functions (fcns) must be equal to number of dofs for setting initial values. %d != %d", len(stg.IniFcn.Fcns), len(stg.IniFcn.Dofs))
	}

	// using file
	if stg.IniFcn.File != "" {
		chk.Panic("initialisation using a file is not implemented yet")
	}

	// loop over functions
	var fcn fun.Func
	for i, fname := range stg.IniFcn.Fcns {

		// get function
		fcn, err = o.Sim.Functions.Get(fname)
		if err != nil {
			return
		}

		// set nodes
		key := stg.IniFcn.Dofs[i]
		for _, nod := range o.Nodes {
			eq := nod.GetEq(key)
			if eq < 0 {
				return chk.Err("dof=%q cannot be found in node=%d for setting initial values", key, nod.Vert.Id)
			}
			o.Sol.Y[eq] = fcn.F(0, nod.Vert.C)
		}
	}
	return
}
