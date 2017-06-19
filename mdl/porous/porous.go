// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package porous implements models for porous media based on the Theory of Porous Media
//  References:
//   [1] Pedroso DM (2015) A consistent u-p formulation for porous media with hysteresis.
//       Int Journal for Numerical Methods in Engineering, 101(8) 606-634
//       http://dx.doi.org/10.1002/nme.4808
//   [2] Pedroso DM (2015) A solution to transient seepage in unsaturated porous media.
//       Computer Methods in Applied Mechanics and Engineering, 285 791-816
//       http://dx.doi.org/10.1016/j.cma.2014.12.009
package porous

import (
	"math"

	"github.com/cpmech/gofem/mdl/conduct"
	"github.com/cpmech/gofem/mdl/fluid"
	"github.com/cpmech/gofem/mdl/retention"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun/dbf"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

// Model holds material parameters for porous media
//  References:
//   [1] Pedroso DM (2015) A consistent u-p formulation for porous media with hysteresis.
//       Int Journal for Numerical Methods in Engineering, 101(8) 606-634
//       http://dx.doi.org/10.1002/nme.4808
//   [2] Pedroso DM (2015) A solution to transient seepage in unsaturated porous media.
//       Computer Methods in Applied Mechanics and Engineering, 285 791-816
//       http://dx.doi.org/10.1016/j.cma.2014.12.009
type Model struct {

	// constants
	NmaxIt  int     // max number iterations in Update
	Itol    float64 // iterations tolerance in Update
	PcZero  float64 // minimum value allowed for pc
	MEtrial bool    // perform Modified-Euler trial to start update process
	ShowR   bool    // show residual values in Update
	AllBE   bool    // use BE for all models, including those that directly implements sl=f(pc)
	Ncns    bool    // use non-consistent method for all derivatives (see [1])
	Ncns2   bool    // use non-consistent method only for second order derivatives (see [1])

	// parameters
	Nf0   float64 // nf0: initial volume fraction of all fluids ~ porosity
	RhoS0 float64 // real (intrinsic) density of solids

	// derived
	Klsat [][]float64 // klsat ÷ Gref
	Kgsat [][]float64 // kgsat ÷ Gref

	// auxiliary models
	Cnd conduct.Model   // liquid-gas conductivity models
	Lrm retention.Model // retention model
	Liq *fluid.Model    // liquid properties
	Gas *fluid.Model    // gas properties

	// auxiliary
	nonrateLrm retention.Nonrate // LRM is of non-rate type
}

// Init initialises this structure
func (o *Model) Init(prms dbf.Params, Cnd conduct.Model, Lrm retention.Model, Liq *fluid.Model, Gas *fluid.Model, grav float64) (err error) {

	// constants
	o.NmaxIt = 20
	o.Itol = 1e-9
	o.PcZero = 1e-10
	o.MEtrial = true

	// read optional constants
	for _, p := range prms {
		switch p.N {
		case "NmaxIt":
			o.NmaxIt = int(p.V)
		case "Itol":
			o.Itol = p.V
		case "PcZero":
			o.PcZero = p.V
		case "MEtrial":
			o.MEtrial = p.V > 0
		case "ShowR":
			o.ShowR = p.V > 0
		case "AllBE":
			o.AllBE = p.V > 0
		case "Ncns":
			o.Ncns = p.V > 0
		case "Ncns2":
			o.Ncns2 = p.V > 0
		}
	}

	// liquid conductivity
	var klx, kly, klz float64
	kl_values, kl_found := prms.GetValues([]string{"klx", "kly", "klz"})
	if !utl.AllTrue(kl_found) {
		p := prms.Find("kl")
		if p == nil {
			return chk.Err("porous model: either 'kl' (isotropic) or ['klx', 'kly', 'klz'] must be given in database of material parameters")
		}
		klx, kly, klz = p.V, p.V, p.V
	} else {
		klx, kly, klz = kl_values[0], kl_values[1], kl_values[2]
	}

	// gas conductivity
	var kgx, kgy, kgz float64
	kg_values, kg_found := prms.GetValues([]string{"kgx", "kgy", "kgz"})
	if !utl.AllTrue(kg_found) {
		p := prms.Find("kg")
		if p == nil {
			return chk.Err("porous model: either 'kg' (isotropic) or ['kgx', 'kgy', 'kgz'] must be given in database of material parameters")
		}
		kgx, kgy, kgz = p.V, p.V, p.V
	} else {
		kgx, kgy, kgz = kg_values[0], kg_values[1], kg_values[2]
	}

	// check conductivities
	KMIN := 1e-14
	if klx < KMIN {
		return chk.Err("porous model: klx (or kl) must be greater than or equal to %g", KMIN)
	}
	if kly < KMIN {
		return chk.Err("porous model: kly (or kl) must be greater than or equal to %g", KMIN)
	}
	if klz < KMIN {
		return chk.Err("porous model: klz (or kl) must be greater than or equal to %g", KMIN)
	}
	if kgx < KMIN {
		return chk.Err("porous modeg: kgx (or kg) must be greater than or equag to %g", KMIN)
	}
	if kgy < KMIN {
		return chk.Err("porous modeg: kgy (or kg) must be greater than or equag to %g", KMIN)
	}
	if kgz < KMIN {
		return chk.Err("porous modeg: kgz (or kg) must be greater than or equag to %g", KMIN)
	}

	// read other paramaters
	prms.Connect(&o.Nf0, "nf0", "porous model")
	prms.Connect(&o.RhoS0, "RhoS0", "porous model")

	// check
	if o.Nf0 < 1e-3 {
		return chk.Err("porous model: porosity nf0 = %g is invalid", o.Nf0)
	}
	if o.RhoS0 < 1e-3 {
		return chk.Err("porous model: intrinsic density of solids RhoS0 = %g is invalid", o.Nf0)
	}
	if grav < 1e-3 {
		return chk.Err("porous model: gravity constant of reference grav = %g is invalid", grav)
	}

	// derived
	o.Klsat = [][]float64{
		{klx / grav, 0, 0},
		{0, kly / grav, 0},
		{0, 0, klz / grav},
	}
	o.Kgsat = [][]float64{
		{kgx / grav, 0, 0},
		{0, kgy / grav, 0},
		{0, 0, kgz / grav},
	}

	// auxiliary models
	if Cnd == nil || Lrm == nil || Liq == nil || Gas == nil {
		return chk.Err("Cnd, Lrm, Liq and Gas models must be all non-nil\n")
	}
	o.Cnd = Cnd
	o.Lrm = Lrm
	o.Liq = Liq
	o.Gas = Gas
	if m, ok := o.Lrm.(retention.Nonrate); ok {
		o.nonrateLrm = m
	}
	return
}

// GetPrms gets (an example) of parameters
func (o Model) GetPrms(example bool) dbf.Params {
	if example {
		return dbf.Params{
			&dbf.P{N: "nf0", V: 0.3},   // [-]
			&dbf.P{N: "RhoS0", V: 2.7}, // [Mg/m³]
			&dbf.P{N: "kl", V: 1e-3},   // [m/s]
			&dbf.P{N: "kg", V: 1e-2},   // [m/s]
		}
	}
	return dbf.Params{
		&dbf.P{N: "nf0", V: o.Nf0},
		&dbf.P{N: "RhoS0", V: o.RhoS0},
	}
}

// NewState creates and initialises a new state structure
//  Note: returns nil on errors
func (o Model) NewState(ρL, ρG, pl, pg float64) (s *State, err error) {
	sl0 := o.Lrm.SlMax()
	sl := sl0
	pc := pg - pl
	if pc > 0 {
		sl, err = retention.Update(o.Lrm, 0, sl0, pc)
		if err != nil {
			return
		}
	}
	ns0 := 1.0 - o.Nf0
	s = &State{ns0, sl, ρL, ρG, 0, false}
	return
}

// Update updates state
//  pl and pg are updated (new) values
func (o Model) Update(s *State, Δpl, Δpg, pl, pg float64) (err error) {

	// auxiliary variables
	slmax := o.Lrm.SlMax()
	slmin := o.Lrm.SlMin()
	Δpc := Δpg - Δpl
	wet := Δpc < 0
	pl0 := pl - Δpl
	pg0 := pg - Δpg
	pc0 := pg0 - pl0
	sl0 := s.A_sl
	pc := pc0 + Δpc
	sl := sl0

	// update liquid saturation
	if pc <= 0.0 {
		sl = slmax // max liquid saturation if capillary pressure is ineffective

	} else if o.nonrateLrm != nil && !o.AllBE {
		sl = o.nonrateLrm.Sl(pc) // handle simple retention models

	} else { // unsaturated case with rate-type model

		// trial liquid saturation update
		fA, e := o.Lrm.Cc(pc0, sl0, wet)
		if e != nil {
			return e
		}
		if o.MEtrial {
			slFE := sl0 + Δpc*fA
			fB, e := o.Lrm.Cc(pc, slFE, wet)
			if e != nil {
				return e
			}
			sl += 0.5 * Δpc * (fA + fB)
		} else {
			sl += Δpc * fA
		}

		// fix trial sl out-of-range values
		if sl < slmin {
			sl = slmin
		}
		if sl > slmax {
			sl = slmax
		}

		// message
		if o.ShowR {
			io.PfYel("%6s%18s%18s%18s%18s%8s\n", "it", "Cc", "sl", "δsl", "r", "ex(r)")
		}

		// backward-Euler update
		var f, r, J, δsl float64
		var it int
		for it = 0; it < o.NmaxIt; it++ {
			f, err = o.Lrm.Cc(pc, sl, wet)
			if err != nil {
				return
			}
			r = sl - sl0 - Δpc*f
			if o.ShowR {
				io.Pfyel("%6d%18.14f%18.14f%18.14f%18.10e%8d\n", it, f, sl, δsl, r, utl.Expon(r))
			}
			if math.Abs(r) < o.Itol {
				break
			}
			J, err = o.Lrm.J(pc, sl, wet)
			if err != nil {
				return
			}
			δsl = -r / (1.0 - Δpc*J)
			sl += δsl
			if math.IsNaN(sl) {
				return chk.Err("NaN found: Δpc=%v f=%v r=%v J=%v sl=%v\n", Δpc, f, r, J, sl)
			}
		}

		// message
		if o.ShowR {
			io.Pfgrey("  pc0=%.6f  sl0=%.6f  Δpl=%.6f  Δpg=%.6f  Δpc=%.6f\n", pc0, sl0, Δpl, Δpg, Δpc)
			io.Pfgrey("  converged with %d iterations\n", it)
		}

		// check convergence
		if it == o.NmaxIt {
			return chk.Err("saturation update failed after %d iterations\n", it)
		}
	}

	// check results
	if pc < 0 && sl < slmax {
		return chk.Err("inconsistent results: saturation must be equal to slmax=%g when the capillary pressure is ineffective. pc = %g < 0 and sl = %g < 1 is incorrect", slmax, pc, sl)
	}
	if sl < slmin {
		return chk.Err("inconsistent results: saturation must be greater than minimum saturation. sl = %g < %g is incorrect", sl, slmin)
	}
	if sl > slmax {
		return chk.Err("inconsistent results: saturation must be smaller than maximum saturation. sl = %g > %g is incorrect", sl, slmax)
	}

	// set state
	s.A_sl = sl             // 2
	s.A_ρL += o.Liq.C * Δpl // 3
	s.A_ρG += o.Gas.C * Δpg // 4
	s.A_Δpc = Δpc           // 5
	s.A_wet = wet           // 6
	return
}

// Ccb (Cc-bar) returns dsl/dpc consistent with the update method
//  See Eq. (54) on page 618 of [1]
func (o Model) Ccb(s *State, pc float64) (dsldpc float64, err error) {
	sl := s.A_sl
	wet := s.A_wet
	Δpc := s.A_Δpc
	f, err := o.Lrm.Cc(pc, sl, wet) // @ n+1
	if err != nil {
		return
	}
	if o.Ncns { // non consistent
		dsldpc = f
		return
	}
	L, err := o.Lrm.L(pc, sl, wet) // @ n+1
	if err != nil {
		return
	}
	J, err := o.Lrm.J(pc, sl, wet) // @ n+1
	if err != nil {
		return
	}
	dsldpc = (f + Δpc*L) / (1.0 - Δpc*J)
	return
}

// Ccd (Cc-dash) returns dCc/dpc consistent with the update method
//  See Eqs. (55) and (56) on page 618 of [1]
func (o Model) Ccd(s *State, pc float64) (dCcdpc float64, err error) {
	sl := s.A_sl
	wet := s.A_wet
	Δpc := s.A_Δpc
	if o.Ncns || o.Ncns2 { // non consistent
		dCcdpc, err = o.Lrm.L(pc, sl, wet) // @ n+1
		return
	}
	f, err := o.Lrm.Cc(pc, sl, wet) // @ n+1
	if err != nil {
		return
	}
	L, Lx, J, Jx, Jy, err := o.Lrm.Derivs(pc, sl, wet)
	if err != nil {
		return
	}
	Ly := Jx
	Ccb := (f + Δpc*L) / (1.0 - Δpc*J)
	LL := Lx + Ly*Ccb
	JJ := Jx + Jy*Ccb
	dCcdpc = (2.0*L + Δpc*LL + (2.0*J+Δpc*JJ)*Ccb) / (1.0 - Δpc*J)
	return
}
