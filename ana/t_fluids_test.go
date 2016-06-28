// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package ana

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func Test_fluids01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("fluids01. properties of water and dry air")

	var water Water
	water.Init()
	io.Pforan("\n>>> water <<<\n")
	io.Pforan("reference temperature: Θ = %23g           [K]\n", water.Θ)
	io.Pforan("bulk modulus @ Θ:      K = %23g           [kPa]\n", water.K)
	io.Pforan("intrinsic density @ Θ: ρ = %23g (0.997)   [Mg/m³]\n", water.Rho)
	io.Pforan("compressibility @ Θ:   C = %23g (4.53e-7) [Mg/(m³・kPa)]\n", water.C)

	var air DryAir
	air.Init()
	io.Pf("\n>>> dry air <<<\n")
	io.Pf("reference temperature: Θ = %23g            [K]\n", air.Θ)
	io.Pf("specific gas constant: R = %23g            [kPa]\n", air.R)
	io.Pf("intrinsic density @ Θ: ρ = %23g (0.001184) [Mg/m³]\n", air.Rho)
	io.Pf("compressibility @ Θ:   C = %23g (1.168e-5) [Mg/(m³・kPa)]\n", air.C)
}
