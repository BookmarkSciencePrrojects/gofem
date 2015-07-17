// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

// Styles
type Styles []plt.Fmt

func GetDefaultStyles(qts Points) Styles {
	sty := make([]plt.Fmt, len(qts))
	for i, q := range qts {
		sty[i].L = io.Sf("x=%v", q.X)
	}
	return sty
}

func GetTexLabel(key, unit string) string {
	l := "$"
	switch key {
	case "time":
		l += "t"
	case "ux":
		l += "u_x"
	case "uy":
		l += "u_y"
	case "uz":
		l += "u_z"
	case "sl":
		l += "s_{\\ell}"
	case "sg":
		l += "s_g"
	case "pl":
		l += "p_{\\ell}"
	case "pg":
		l += "p_g"
	case "sx":
		l += "\\sigma_x"
	case "sy":
		l += "\\sigma_y"
	case "sz":
		l += "\\sigma_z"
	case "sxy":
		l += "\\sigma_{xy}"
	case "syz":
		l += "\\sigma_{yz}"
	case "szx":
		l += "\\sigma_{zx}"
	case "ex_nwlx", "nwlx":
		l += "n_{\\ell}\\cdot w_{\\ell x}"
	case "ex_nwly", "nwly":
		l += "n_{\\ell}\\cdot w_{\\ell y}"
	case "ex_nwlz", "nwlz":
		l += "n_{\\ell}\\cdot w_{\\ell z}"
	case "ompb":
		l += "\\bar{\\omega}_p"
	default:
		l += key
	}
	if unit != "" {
		l += "\\;" + unit
	}
	l += "$"
	return l
}
