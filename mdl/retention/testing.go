// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package retention

import (
	"math"
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

// checkDerivs checks derivatives
func checkDerivs(tst *testing.T, mdl Model, pc0, sl0, pcf float64, npts int, tolCc, tolD1a, tolD1b, tolD2a, tolD2b float64, verbose bool, pcSkip []float64, tolSkip float64, doplot bool) {

	// nonrate model
	nr_mdl, is_nonrate := mdl.(Nonrate)
	io.Pforan("is_nonrate = %v\n", is_nonrate)

	// for all pc stations
	Pc := utl.LinSpace(pc0, pcf, npts)
	Sl := make([]float64, npts)
	Sl[0] = sl0
	var err error
	for i := 1; i < npts; i++ {

		// update and plot
		Sl[i], err = Update(mdl, Pc[i-1], Sl[i-1], Pc[i]-Pc[i-1])
		if err != nil {
			tst.Errorf("Update failed: %v\n", err)
			return
		}
		//if doplot {
		//plt.PlotOne(Pc[i], Sl[i], "'ko', clip_on=0")
		//}

		// skip point on checking of derivatives
		if doskip(Pc[i], pcSkip, tolSkip) {
			continue
		}

		// wetting flag
		wet := Pc[i]-Pc[i-1] < 0

		// check Cc = dsl/dpc
		io.Pforan("\npc=%g, sl=%g, wetting=%v\n", Pc[i], Sl[i], wet)
		if is_nonrate {

			// analytical Cc
			Cc_ana, err := mdl.Cc(Pc[i], Sl[i], wet)
			if err != nil {
				tst.Errorf("Cc failed: %v\n", err)
				return
			}

			// numerical Cc
			chk.DerivScaSca(tst, "Cc = ∂sl/∂pc    ", tolCc, Cc_ana, Pc[i], 1e-3, verbose, func(x float64) (float64, error) {
				return nr_mdl.Sl(x), nil
			})
		}

		// compute all derivatives
		L, Lx, J, Jx, Jy, err := mdl.Derivs(Pc[i], Sl[i], wet)
		if err != nil {
			tst.Errorf("Derivs failed: %v\n", err)
			return
		}
		L_ana_A := L
		L_ana_B, err := mdl.L(Pc[i], Sl[i], wet)
		if err != nil {
			tst.Errorf("L failed: %v\n", err)
			return
		}
		Lx_ana := Lx
		Jx_ana := Jx
		Jy_ana := Jy
		J_ana_A := J
		J_ana_B, err := mdl.J(Pc[i], Sl[i], wet)
		if err != nil {
			tst.Errorf("J failed: %v\n", err)
			return
		}

		// numerical L = ∂Cc/∂pc
		chk.DerivScaSca(tst, "L  = ∂Cc/∂pc    ", tolD1a, L_ana_A, Pc[i], 1e-3, verbose, func(x float64) (float64, error) {
			Cctmp, _ := mdl.Cc(x, Sl[i], wet)
			return Cctmp, nil
		})

		// numerical Lx := ∂²Cc/∂pc²
		chk.DerivScaSca(tst, "Lx = ∂²Cc/∂pc²  ", tolD2a, Lx_ana, Pc[i], 1e-3, verbose, func(x float64) (float64, error) {
			Ltmp, _, _, _, _, _ := mdl.Derivs(x, Sl[i], wet)
			return Ltmp, nil
		})

		// numerical J := ∂Cc/∂sl (version A)
		chk.DerivScaSca(tst, "J  = ∂Cc/∂sl    ", tolD1b, J_ana_A, Sl[i], 1e-3, verbose, func(x float64) (float64, error) {
			Ccval, _ := mdl.Cc(Pc[i], x, wet)
			return Ccval, nil
		})

		// numerical Jx := ∂²Cc/(∂pc ∂sl)
		chk.DerivScaSca(tst, "Jx = ∂²Cc/∂pc∂sl", tolD2b, Jx_ana, Sl[i], 1e-3, verbose, func(x float64) (float64, error) {
			Ltmp, _, _, _, _, _ := mdl.Derivs(Pc[i], x, wet)
			return Ltmp, nil
		})

		// numerical Jy := ∂²Cc/∂sl²
		chk.DerivScaSca(tst, "Jy = ∂²Cc/∂sl²  ", tolD2b, Jy_ana, Sl[i], 1e-3, verbose, func(x float64) (float64, error) {
			Jtmp, _ := mdl.J(Pc[i], x, wet)
			return Jtmp, nil
		})

		// check A and B derivatives
		chk.Float64(tst, "L_A == L_B", 1e-17, L_ana_A, L_ana_B)
		chk.Float64(tst, "J_A == J_B", 1e-17, J_ana_A, J_ana_B)
	}
}

// doskip analyse whether a point should be skip or not
func doskip(x float64, xskip []float64, tol float64) bool {
	for _, v := range xskip {
		if math.Abs(x-v) < tol {
			return true
		}
	}
	return false
}
