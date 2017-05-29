// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package porous

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

// Driver run simulations with models for porous media
type Driver struct {

	// input
	Mdl *Model // porous model

	// settings
	Silent bool    // do not show error messages
	TolCcb float64 // tolerance to check Ccb
	TolCcd float64 // tolerance to check Ccd
	VerD   bool    // verbose check of D

	// check D matrix
	TstD *testing.T // if != nil, do check consistent matrix

	// results
	Res []*State // results
}

// Init initialises driver
func (o *Driver) Init(mdl *Model) (err error) {
	o.Mdl = mdl
	o.TolCcb = 1e-7
	o.TolCcd = 1e-7
	o.VerD = chk.Verbose
	return
}

// Run runs simulation
func (o *Driver) Run(Pc []float64) (err error) {

	// allocate results arrays
	np := len(Pc)
	o.Res = make([]*State, np)

	// initialise first state
	o.Res[0], err = o.Mdl.NewState(o.Mdl.Liq.R0, o.Mdl.Gas.R0, -Pc[0], 0)
	if err != nil {
		return
	}

	// update states
	var pcOld, pcNew, Δpc, Ccb, Ccbtmp, Ccd float64
	var stmp State
	for i := 1; i < np; i++ {

		// increment
		pcOld = Pc[i-1]
		pcNew = Pc[i]
		Δpc = pcNew - pcOld

		// update
		o.Res[i] = o.Res[i-1].GetCopy()
		err = o.Mdl.Update(o.Res[i], -Δpc, 0, -pcNew, 0)
		if err != nil {
			return
		}

		// check consistent moduli
		if o.TstD != nil {

			// check Ccb
			Ccb, err = o.Mdl.Ccb(o.Res[i], pcNew)
			if err != nil {
				return
			}
			chk.DerivScaSca(o.TstD, io.Sf("Ccb @ %.3f,%.4f", pcNew, o.Res[i].A_sl), o.TolCcb, Ccb, pcNew, 1e-3, o.VerD, func(x float64) (float64, error) {
				Δpc = x - pcOld
				stmp.Set(o.Res[i-1])
				e := o.Mdl.Update(&stmp, -Δpc, 0, -x, 0)
				return stmp.A_sl, e
			})

			// check Ccd
			Ccd, err = o.Mdl.Ccd(o.Res[i], pcNew)
			if err != nil {
				return
			}
			chk.DerivScaSca(o.TstD, io.Sf("Ccd @ %.3f,%.4f", pcNew, o.Res[i].A_sl), o.TolCcd, Ccd, pcNew, 1e-3, o.VerD, func(x float64) (float64, error) {
				Δpc = x - pcOld
				stmp.Set(o.Res[i-1])
				e := o.Mdl.Update(&stmp, -Δpc, 0, -x, 0)
				if e != nil {
					return 0, e
				}
				Ccbtmp, e = o.Mdl.Ccb(&stmp, x)
				return Ccbtmp, e
			})
		}
	}
	return
}
