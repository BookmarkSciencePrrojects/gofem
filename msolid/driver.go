// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/num"
)

// Driver run simulations with constitutive models for solids
type Driver struct {

	// input
	nsig  int   // number of stress components
	model Model // solid model

	// settings
	Silent  bool    // do not show error messages
	CheckD  bool    // do check consistent matrix
	UseDfwd bool    // use DerivFwd (forward differences) instead of DerivCen (central differences) when checking D
	TolD    float64 // tolerance to check consistent matrix
	VerD    bool    // verbose check of D
	WithPC  bool    // with predictor-corrector data

	// results
	Res []*State    // stress/ivs results
	Eps [][]float64 // strains

	// for checking consistent matrix
	D [][]float64 // consistent matrix

	// for predictor-corrector plots
	PreCor [][]float64 // predictor-corrector stresses
}

// Init initialises driver
func (o *Driver) Init(simfnk, modelname string, ndim int, pstress bool, prms fun.Prms) (err error) {
	o.nsig = 2 * ndim
	getnew := false
	o.model = GetModel(simfnk, "solidmat", modelname, getnew)
	if o.model == nil {
		return chk.Err("cannot get model named %q\n", modelname)
	}
	err = o.model.Init(ndim, pstress, prms)
	if err != nil {
		return
	}
	o.D = la.MatAlloc(o.nsig, o.nsig)
	o.TolD = 1e-8
	o.VerD = chk.Verbose
	return
}

// Run runs simulation
func (o *Driver) Run(pth *Path) (err error) {

	// specialised models
	var sml Small
	var eup SmallStrainUpdater
	switch m := o.model.(type) {
	case Small:
		sml = m
	case SmallStrainUpdater:
		eup = m
	default:
		return chk.Err("cannot handle large-deformation models yet\n")
	}

	// elastoplastic model
	epm := o.model.(EPmodel)

	// initial stresses
	σ := make([]float64, o.nsig)
	σ[0] = pth.MultS * pth.Sx[0]
	σ[1] = pth.MultS * pth.Sy[0]
	σ[2] = pth.MultS * pth.Sz[0]

	// allocate results arrays
	nr := 1 + (pth.Size()-1)*pth.Nincs
	if nr < 2 {
		return chk.Err(_driver_err04, pth.Size(), pth.Nincs)
	}
	o.Res = make([]*State, nr)
	o.Eps = la.MatAlloc(nr, o.nsig)
	for i := 0; i < nr; i++ {
		o.Res[i], err = o.model.InitIntVars(σ)
		if err != nil {
			return
		}
	}

	// auxiliary variables
	Δσ := make([]float64, o.nsig)
	Δε := make([]float64, o.nsig)

	// variables for checking D
	var tmp float64
	var εold, εnew, Δεtmp []float64
	var stmp *State
	derivfcn := num.DerivCen
	if o.CheckD {
		εold = make([]float64, o.nsig)
		εnew = make([]float64, o.nsig)
		Δεtmp = make([]float64, o.nsig)
		stmp, err = o.model.InitIntVars(σ)
		if err != nil {
			return
		}
		if o.UseDfwd {
			derivfcn = num.DerivFwd
		}
	}

	// update states
	k := 1
	for i := 1; i < pth.Size(); i++ {

		// stress path
		if pth.UseS[i] > 0 {
			Δσ[0] = pth.MultS * (pth.Sx[i] - pth.Sx[i-1]) / float64(pth.Nincs)
			Δσ[1] = pth.MultS * (pth.Sy[i] - pth.Sy[i-1]) / float64(pth.Nincs)
			Δσ[2] = pth.MultS * (pth.Sz[i] - pth.Sz[i-1]) / float64(pth.Nincs)
			for inc := 0; inc < pth.Nincs; inc++ {

				// update
				o.Res[k].Set(o.Res[k-1])
				copy(o.Eps[k], o.Eps[k-1])
				if eup != nil {
					err = eup.StrainUpdate(o.Res[k], Δσ)
				}
				if err != nil {
					if !o.Silent {
						io.Pfred(_driver_err01, err)
					}
					return
				}
				k += 1
			}
		}

		// strain path
		if pth.UseE[i] > 0 {
			Δε[0] = pth.MultE * (pth.Ex[i] - pth.Ex[i-1]) / float64(pth.Nincs)
			Δε[1] = pth.MultE * (pth.Ey[i] - pth.Ey[i-1]) / float64(pth.Nincs)
			Δε[2] = pth.MultE * (pth.Ez[i] - pth.Ez[i-1]) / float64(pth.Nincs)
			for inc := 0; inc < pth.Nincs; inc++ {

				// update strains
				la.VecAdd2(o.Eps[k], 1, o.Eps[k-1], 1, Δε) // εnew = εold + Δε

				// update stresses
				o.Res[k].Set(o.Res[k-1])
				err = sml.Update(o.Res[k], o.Eps[k], Δε)
				if err != nil {
					if !o.Silent {
						io.Pfred(_driver_err02, err)
					}
					return
				}
				if epm != nil {
					tmp := o.Res[k-1].GetCopy()
					epm.ElastUpdate(tmp, Δε)
					o.PreCor = append(o.PreCor, tmp.Sig, o.Res[k].Sig)
				}

				// check consistent matrix
				if o.CheckD {
					firstIt := false
					err = sml.CalcD(o.D, o.Res[k], firstIt)
					if err != nil {
						return chk.Err(_driver_err03, err)
					}
					copy(εold, o.Eps[k-1])
					copy(εnew, o.Eps[k])
					has_error := false
					if o.VerD {
						io.Pf("\n")
					}
					for i := 0; i < o.nsig; i++ {
						for j := 0; j < o.nsig; j++ {
							dnum := derivfcn(func(x float64, args ...interface{}) (res float64) {
								tmp, εnew[j] = εnew[j], x
								for l := 0; l < o.nsig; l++ {
									Δεtmp[l] = εnew[l] - εold[l]
								}
								stmp.Set(o.Res[k-1])
								err = sml.Update(stmp, εnew, Δεtmp)
								res, εnew[j] = stmp.Sig[i], tmp
								return
							}, εnew[j])
							err := chk.PrintAnaNum(io.Sf("D[%d][%d]", i, j), o.TolD, o.D[i][j], dnum, o.VerD)
							if err != nil {
								has_error = true
							}
						}
					}
					if has_error {
						return chk.Err(_driver_err03, "ana-num comparison failed\n")
					}
				}
				k += 1
			}
		}
	}
	return
}

// error messages
var (
	_driver_err01 = "strain update failed\n%v\n"
	_driver_err02 = "stress update failed\n%v\n"
	_driver_err03 = "check of consistent matrix failed:\n %v\n\n"
	_driver_err04 = "size of path is incorrect. Size=%d, Nincs=%d\n"
)
