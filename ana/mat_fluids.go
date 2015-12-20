// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package ana

// Water handles the properties of water
type Water struct {
	Θ   float64 // reference temperature; default = 25°C or 298.15K
	K   float64 // bulk modulus @ reference temperature
	Rho float64 // intrinsic density @ reference temperature
	C   float64 // compressibility @ reference temperature
}

// DryAir handles the properties of dry air
type DryAir struct {
	Θ    float64 // reference temperature; default = 25°C or 298.15K
	R    float64 // specific ideal gas constant
	Patm float64 // absolute atmospheric pressure
	Rho  float64 // intrinsic density @ reference temperature
	C    float64 // compressibility @ reference temperature
}

// Init initialises data
func (o *Water) Init() {
	o.Θ = 298.15      // [K]      25°C
	o.K = 2.2e6       // [kPa]    25°C
	o.Rho = 0.9970479 // [Mg/m³]  25°C
	o.C = o.Rho / o.K // [Mg/(m³・kPa)]
}

// Init initialises data
func (o *DryAir) Init() {
	o.Θ = 298.15                 // [K]          25°C
	o.R = 287.058                // [J/(kg・K)]  25°C  [kJ/(Mg・K)]  [kPa・m³/(Mg・K)]
	o.Patm = 101.325             // [kPa]
	o.Rho = o.Patm / (o.R * o.Θ) // [Mg/m³]      25°C
	o.C = 1.0 / (o.R * o.Θ)      // [Mg/(m³・kPa)]
}
