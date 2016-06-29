// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package seepage

import (
	"math"

	"github.com/cpmech/gosl/io"
)

func GetSeepFaceFlags(extra string) (Macaulay bool, BetRamp, Kappa float64) {

	// defaults
	Macaulay = false
	BetRamp = math.Ln2 / 0.01
	Kappa = 1.0

	// use macaulay function ?
	if s_mac, found := io.Keycode(extra, "mac"); found {
		Macaulay = io.Atob(s_mac)
	}

	// coefficient for smooth ramp function
	if s_bet, found := io.Keycode(extra, "bet"); found {
		BetRamp = io.Atof(s_bet)
	}

	// κ coefficient
	if s_kap, found := io.Keycode(extra, "kap"); found {
		Kappa = io.Atof(s_kap)
	}
	return
}
