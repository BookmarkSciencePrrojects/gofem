// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package ana

import (
	"math"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

// CrossSection computes cross-sectional moments of inertia and other properties
//
//   typ : rectangle
//         circle                             tw
//         I-beam                         -->| |<--
//                                    ___    | |     ___
//   ^ s       +-------+            tf |   ########   |
//   |         |       |              ---  ########   |
//   |         |       |                      ##      |
//   +----> r  |       | h = hei              ##      | h = hei
//             |       |                      ##      |
//             |       |              ---  ########   |
//             +-------+            tf_|_  ########  ---
//              b = wid                    b = wid
//
type CrossSection struct {

	// input
	Type string  // "rectangle", "I-beam" or "circle"
	Unit string  // unit of length
	Wid  float64 // width (b) if not circular
	Hei  float64 // height (h) if not circular
	Tf   float64 // flange thickness if I-beam
	Tw   float64 // web thickness if I-beam
	R    float64 // radius if circular

	// derived
	A   float64 // cross-sectional area
	Irr float64 // major cross-section moment of inertia (about r-axis)
	Iss float64 // minor cross-section moment of inertia (about s-axis)
	Jtt float64 // torsional constant
}

// Init initialises structure and computes moment of inertia
func (o *CrossSection) Init(typ, unitLen string, wid, hei, tf, tw, rad float64) {

	// input data
	o.Type, o.Unit, o.Wid, o.Hei, o.Tf, o.Tw, o.R = typ, unitLen, wid, hei, tf, tw, rad

	// derived
	switch typ {
	case "rectangle":
		b, h := wid, hei
		b3 := b * b * b
		h3 := h * h * h
		o.A = b * h
		o.Irr = b * h3 / 12.0
		o.Iss = b3 * h / 12.0
		if b == h {
			o.Jtt = 9.0 * b3 * b / 64.0
		} else {
			if b > h {
				b, h = h, b
			}
			o.Jtt = h * b3 * (1.0/3.0 - 0.21*(b/h)*(1.0-b*b3/(12.0*h*h3))) // approximate
		}

	case "I-beam":
		b, h := wid, hei
		b3 := b * b * b
		h3 := h * h * h
		tf3 := tf * tf * tf
		tw3 := tw * tw * tw
		l := h - 2.0*tf
		l3 := l * l * l
		o.A = b*h - l*(b-tw)
		o.Irr = b*h3/12.0 - (b-tw)*l3/12.0
		o.Iss = l*tw3/12.0 + tf*b3/6.0
		o.Jtt = (2.0*b*tf3 + (h-2.0*tf)*tw3) / 3.0

	case "circle":
		r2 := rad * rad
		o.A = math.Pi * r2
		o.Irr = math.Pi * r2 * r2 / 4.0
		o.Iss = o.Irr
		o.Jtt = o.Irr + o.Iss

	default:
		chk.Panic("cross-section type %q is unavailable", typ)
	}
}

// GetMatString returns string representation of cross-section for .mat file
func (o *CrossSection) GetMatString(numfmt string) string {
	l := io.Sf("        {\"n\":\"A\",   \"v\":"+numfmt+", \"u\":\""+o.Unit+"²\"},\n", o.A)
	l += io.Sf("        {\"n\":\"Irr\", \"v\":"+numfmt+", \"u\":\""+o.Unit+"⁴\"},\n", o.Irr)
	l += io.Sf("        {\"n\":\"Iss\", \"v\":"+numfmt+", \"u\":\""+o.Unit+"⁴\"},\n", o.Iss)
	l += io.Sf("        {\"n\":\"Jtt\", \"v\":"+numfmt+", \"u\":\""+o.Unit+"⁴\"}", o.Jtt)
	return l
}

// Material holds parameters of some reference materials
type Material struct {

	// input
	Type     string // type of material; e.g. "steel"
	UnitPres string // unit of pressure

	// derived
	UnitDens string  // unit of density
	Desc     string  // description
	E        float64 // Young's modulus
	Nu       float64 // Poisson's coefficient
	G        float64 // shear modulus
	Rho      float64 // density
}

// Init initialises material paramters
//  Input:
//   unitPres:  "kPa" => E:[kPa], rho:[Mg/m^3]
//  		    "MPa" => E:[MPa], rho:[Gg/m^3]
//  		    "GPa" => E:[GPa], rho:[Tg/m^3]
func (o *Material) Init(typ, unitPres string) {

	// material data
	switch typ {
	case "steel":
		o.Desc = "Steel: structural A36"
		o.E = 200000.0  // [MPa]
		o.Nu = 0.32     // [-]
		o.Rho = 7.85e-3 // [Gg/m³]
	case "aluminum":
		o.Desc = "Aluminum: 2014-T6"
		o.E = 73100.0   // [MPa]
		o.Nu = 0.35     // [-]
		o.Rho = 2.79e-3 // [Gg/m³]
	case "concrete-low":
		o.Desc = "Concrete: low strength"
		o.E = 22100.0   // [MPa]
		o.Nu = 0.15     // [-]
		o.Rho = 2.38e-3 // [Gg/m³]
	case "concrete-high":
		o.Desc = "Concrete: high strength"
		o.E = 30000.0   // [MPa]
		o.Nu = 0.15     // [-]
		o.Rho = 2.38e-3 // [Gg/m³]
	case "soft-soil":
		o.Desc = "Soil: soft"
		o.E = 10.0      // [MPa]
		o.Nu = 0.30     // [-]
		o.Rho = 1.80e-3 // [Gg/m³]
	case "wood-douglas-fir":
		o.Desc = "Wood: Douglas-fir"
		o.E = 13100.0   // [MPa]
		o.Nu = 0.29     // [-]
		o.Rho = 4.70e-4 // [Gg/m³]
	default:
		chk.Panic("material type %q is unavailable", typ)
	}

	// set unit
	o.UnitPres = unitPres
	MPa_to_unitPres := 1.0   // convert from MPa to unitPress (e.g. kPa)
	GgByM3_toUnitDens := 1.0 // convert from Gg/m3 to unitPress (e.g. Mg/m³)
	switch unitPres {
	case "kPa":
		o.UnitDens = "Mg/m³"
		MPa_to_unitPres = 1e3   // convert from MPa to kPa
		GgByM3_toUnitDens = 1e3 // convert from Gg/m3 to Mg/m³
	case "MPa":
		o.UnitDens = "Gg/m³"
	case "GPa":
		o.UnitDens = "Tg/m³"
		MPa_to_unitPres = 1e-3   // convert from MPa to GPa
		GgByM3_toUnitDens = 1e-3 // convert from Gg/m3 to Tg/m³
	default:
		chk.Panic("unit of pressure %q is invalid", unitPres)
	}

	// convert values to requested units
	o.E = o.E * MPa_to_unitPres
	o.Rho = o.Rho * GgByM3_toUnitDens

	// derived quantity
	o.G = o.E / (2.0 * (1.0 + o.Nu))
}

func (o *Material) GetMatString(name, model, numfmt string, section *CrossSection) string {
	l := io.Sf("      \"name\" : %q,\n", name)
	if model != "" {
		l += io.Sf("      \"model\" : %q,\n", model)
	}
	l += "      \"prms\" : [\n"
	l += io.Sf("        {\"n\":\"E\",   \"v\":"+numfmt+", \"u\":%q},\n", o.E, o.UnitPres)
	l += io.Sf("        {\"n\":\"G\",   \"v\":"+numfmt+", \"u\":%q},\n", o.G, o.UnitPres)
	l += io.Sf("        {\"n\":\"nu\",  \"v\":"+numfmt+", \"u\":%q},\n", o.Nu, "-")
	l += io.Sf("        {\"n\":\"rho\", \"v\":"+numfmt+", \"u\":%q}", o.Rho, o.UnitDens)
	if section != nil {
		l += ",\n" + section.GetMatString(numfmt)
	}
	l += io.Sf("\n      ]")
	return l
}
