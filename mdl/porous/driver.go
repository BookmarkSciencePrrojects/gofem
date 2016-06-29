// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package porous

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/num"
)

// Driver run simulations with models for porous media
type Driver struct {

	// input
	Mdl *Model // porous model

	// settings
	Silent  bool    // do not show error messages
	CheckD  bool    // do check consistent matrix
	UseDfwd bool    // use DerivFwd (forward differences) instead of DerivCen (central differences) when checking D
	TolCcb  float64 // tolerance to check Ccb
	TolCcd  float64 // tolerance to check Ccd
	VerD    bool    // verbose check of D

	// results
	Res []*State // results
}

// Init initialises driver
func (o *Driver) Init(mdl *Model) (err error) {
	o.Mdl = mdl
	o.TolCcb = 1e-7
	o.TolCcd = 1e-7
	o.VerD = chk.Verbose
	o.CheckD = true
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

	// auxiliary
	derivfcn := num.DerivCen
	if o.UseDfwd {
		derivfcn = num.DerivFwd
	}

	// update states
	var pcOld, pcNew, Δpc, tmp, Ccb, Ccbtmp, Ccd float64
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
		if o.CheckD {

			// check Ccb
			Ccb, err = o.Mdl.Ccb(o.Res[i], pcNew)
			if err != nil {
				return
			}
			var has_errors bool
			dnum := derivfcn(func(x float64, args ...interface{}) (res float64) {
				tmp, pcNew = pcNew, x
				Δpc = pcNew - pcOld
				stmp.Set(o.Res[i-1])
				e := o.Mdl.Update(&stmp, -Δpc, 0, -pcNew, 0)
				if e != nil {
					has_errors = true
				}
				res, pcNew = stmp.A_sl, tmp
				return
			}, pcNew)
			if has_errors {
				return chk.Err("problems arised during update in numerical derivative for Ccb")
			}
			err = chk.PrintAnaNum(io.Sf("Ccb @ %.3f,%.4f", pcNew, o.Res[i].A_sl), o.TolCcb, Ccb, dnum, o.VerD)
			if err != nil {
				return
			}

			// check Ccd
			Ccd, err = o.Mdl.Ccd(o.Res[i], pcNew)
			if err != nil {
				return
			}
			has_errors = false
			dnum = derivfcn(func(x float64, args ...interface{}) (res float64) {
				tmp, pcNew = pcNew, x
				Δpc = pcNew - pcOld
				stmp.Set(o.Res[i-1])
				e := o.Mdl.Update(&stmp, -Δpc, 0, -pcNew, 0)
				if e != nil {
					has_errors = true
				}
				Ccbtmp, _ = o.Mdl.Ccb(&stmp, pcNew)
				res, pcNew = Ccbtmp, tmp
				return
			}, pcNew)
			if has_errors {
				return chk.Err("problems arised during update in numerical derivative for Ccd")
			}
			err = chk.PrintAnaNum(io.Sf("Ccd @ %.3f,%.4f", pcNew, o.Res[i].A_sl), o.TolCcd, Ccd, dnum, o.VerD)
			if err != nil {
				return
			}
		}
	}
	return
}
