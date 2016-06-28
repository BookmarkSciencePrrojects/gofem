// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package por

// State holds state variables for porous media with liquid and gas
//  References:
//   [1] Pedroso DM (2015) A consistent u-p formulation for porous media with hysteresis.
//       Int Journal for Numerical Methods in Engineering, 101(8) 606-634
//       http://dx.doi.org/10.1002/nme.4808
//   [2] Pedroso DM (2015) A solution to transient seepage in unsaturated porous media.
//       Computer Methods in Applied Mechanics and Engineering, 285 791-816
//       http://dx.doi.org/10.1016/j.cma.2014.12.009
type State struct {
	A_ns0 float64 // 1 initial partial fraction of solids
	A_sl  float64 // 2 liquid saturation
	A_ρL  float64 // 3 real (intrinsic) density of liquid
	A_ρG  float64 // 4 real (intrinsic) density of gas
	A_Δpc float64 // 5 step increment of capillary pressure
	A_wet bool    // 6 wetting flag
}

// GetCopy returns a copy of State
func (o State) GetCopy() *State {
	return &State{
		o.A_ns0, // 1
		o.A_sl,  // 2
		o.A_ρL,  // 3
		o.A_ρG,  // 4
		o.A_Δpc, // 5
		o.A_wet, // 6
	}
}

// Set sets this State with another State
func (o *State) Set(s *State) {
	o.A_ns0 = s.A_ns0 // 1
	o.A_sl = s.A_sl   // 2
	o.A_ρL = s.A_ρL   // 3
	o.A_ρG = s.A_ρG   // 4
	o.A_Δpc = s.A_Δpc // 5
	o.A_wet = s.A_wet // 6
}

// LsVars hold data for liquid-solid computations
type LsVars struct {
	A_ρl, A_ρ, A_p, Cpl, Cvs                float64
	Dρdpl, Dpdpl, DCpldpl, DCvsdpl, Dklrdpl float64
	DρldusM, DρdusM, DCpldusM               float64
}

// LgsVars hold data for liquid-gas-solid computations
type LgsVars struct {
	A_ρl, A_ρg       float64
	A_ρ, A_p         float64
	Cpl, Cpg, Cvs    float64
	Dpl, Dpg, Dvs    float64
	Dklrdpl, Dklrdpg float64
	Dkgrdpl, Dkgrdpg float64
	DρlduM, DρgduM   float64
	Dρdpl, Dρdpg     float64
	DρduM            float64
	Dpdpl, Dpdpg     float64
	DCpldpl, DCpldpg float64
	DCpgdpl, DCpgdpg float64
	DCvsdpl, DCvsdpg float64
	DDpldpl, DDpldpg float64
	DDpgdpl, DDpgdpg float64
	DDvsdpl, DDvsdpg float64
	DCplduM, DCpgduM float64
	DDplduM, DDpgduM float64
}

// CalcLs calculates variables for liquid-solid simulations
func (o Model) CalcLs(res *LsVars, sta *State, pl, divu float64, derivs bool) (err error) {

	// auxiliary
	ns0 := sta.A_ns0
	sl := sta.A_sl
	ρL := sta.A_ρL
	Cl := o.Liq.C
	ρS := o.RhoS0

	// n variables; Eqs (13) and (28) of [1]
	ns := (1.0 - divu) * ns0
	nf := 1.0 - ns
	nl := nf * sl

	// ρ variables; Eq (13) of [1]
	ρs := ns * ρS
	res.A_ρl = nl * ρL
	res.A_ρ = res.A_ρl + ρs

	// capillary pressure and pore-fluid pressure
	pc := -pl
	res.A_p = pl * sl // Eq. (16) of [1]

	// moduli
	Ccb, e := o.Ccb(sta, pc)
	if e != nil {
		return e
	}
	res.Cpl = nf * (sl*Cl - ρL*Ccb) // Eq (32a) of [1]
	res.Cvs = sl * ρL               // Eq (32b) of [1]

	// derivatives
	if derivs {

		// Ccd
		Ccd, e := o.Ccd(sta, pc)
		if e != nil {
			return e
		}

		// derivatives w.r.t pl
		res.Dρdpl = nf * (sl*Cl - ρL*Ccb)        // Eq (A.9) of [1]
		res.Dpdpl = sl + pc*Ccb                  // Eq (A.11) of [1]
		res.DCpldpl = nf * (ρL*Ccd - 2.0*Ccb*Cl) // Eq (A.2) of[1]
		res.DCvsdpl = sl*Cl - Ccb*ρL             // Eq (A.4) of [1]
		res.Dklrdpl = -o.Cnd.DklrDsl(sl) * Ccb   // Eq (A.7) of [1]

		// derivatives w.r.t us (multipliers only)
		res.DρldusM = sl * ρL * ns0
		res.DρdusM = (sl*ρL - ρS) * ns0       // Eq (A.10) of [1]
		res.DCpldusM = (sl*Cl - ρL*Ccb) * ns0 // Eq (A.3) of [1]
	}
	return
}

// CalcLgs calculates variables for liquid-gas-solid simulations
func (o Model) CalcLgs(res *LgsVars, sta *State, pl, pg, divus float64, derivs bool) (err error) {

	// auxiliary
	ns0 := sta.A_ns0
	sl := sta.A_sl
	sg := 1.0 - sl
	ρL := sta.A_ρL
	ρG := sta.A_ρG
	Cl := o.Liq.C
	Cg := o.Gas.C
	ρS := o.RhoS0

	// n variables
	ns := (1.0 - divus) * ns0
	nf := 1.0 - ns
	nl := nf * sl
	ng := nf * sg

	// ρ variables
	ρs := ns * ρS
	res.A_ρl = nl * ρL
	res.A_ρg = ng * ρG
	res.A_ρ = res.A_ρl + res.A_ρg + ρs

	// capillary pressure and pore-fluid pressure
	pc := pg - pl
	res.A_p = pl*sl + pg*sg

	// moduli
	Cc, e := o.Ccb(sta, pc)
	if e != nil {
		return e
	}
	res.Cpl = nf * (sl*Cl - ρL*Cc)
	res.Cpg = nf * ρL * Cc
	res.Cvs = sl * ρL
	res.Dpl = nf * ρG * Cc
	res.Dpg = nf * (sg*Cg - ρG*Cc)
	res.Dvs = sg * ρG

	// derivatives
	if derivs {

		// Ccd
		Ccd, e := o.Ccd(sta, pc)
		if e != nil {
			return e
		}

		// conductivity multipliers
		dklrdsl := o.Cnd.DklrDsl(sl)
		dkgrdsg := o.Cnd.DkgrDsg(sg)
		res.Dklrdpl = -dklrdsl * Cc
		res.Dklrdpg = dklrdsl * Cc
		res.Dkgrdpl = dkgrdsg * Cc
		res.Dkgrdpg = -dkgrdsg * Cc

		// partial densities
		res.DρlduM = sl * ρL * ns0
		res.DρgduM = sg * ρG * ns0

		// mixture density
		res.Dρdpl = nf * (sl*Cl - ρL*Cc + ρG*Cc)
		res.Dρdpg = nf * (sg*Cg - ρG*Cc + ρL*Cc)
		res.DρduM = (sl*ρL + sg*ρG - ρS) * ns0

		// pressure in pores
		res.Dpdpl = sl + pc*Cc
		res.Dpdpg = sg - pc*Cc

		// derivatives of C coefficients
		res.DCpldpl = nf * (ρL*Ccd - 2.0*Cc*Cl)
		res.DCpldpg = nf * (Cc*Cl - ρL*Ccd)
		res.DCpgdpl = nf * (Cl*Cc - ρL*Ccd)
		res.DCpgdpg = nf * ρL * Ccd
		res.DCvsdpl = sl*Cl - Cc*ρL
		res.DCvsdpg = Cc * ρL

		// derivatives of D coefficients
		res.DDpldpl = -nf * ρG * Ccd
		res.DDpldpg = nf * (ρG*Ccd + Cg*Cc)
		res.DDpgdpl = nf * (Cc*Cg + ρG*Ccd)
		res.DDpgdpg = -nf * (ρG*Ccd + 2.0*Cg*Cc)
		res.DDvsdpl = Cc * ρG
		res.DDvsdpg = sg*Cg - Cc*ρG

		// derivatives w.r.t u (multiplier)
		res.DCplduM = (sl*Cl - ρL*Cc) * ns0
		res.DCpgduM = ρL * Cc * ns0
		res.DDplduM = ρG * Cc * ns0
		res.DDpgduM = (sg*Cg - ρG*Cc) * ns0
	}
	return
}
