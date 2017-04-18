// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

// Styles
type Styles []plt.A

func GetDefaultStyles(qts Points) Styles {
	sty := make([]plt.A, len(qts))
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
	case "pc":
		l += "p_c"
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
	case "ex_nwgx", "nwgx":
		l += "n_{g}\\cdot w_{g x}"
	case "ex_nwgy", "nwgy":
		l += "n_{g}\\cdot w_{g y}"
	case "ex_nwgz", "nwgz":
		l += "n_{g}\\cdot w_{g z}"
	case "RhoL":
		l += "\\rho^{\\ell}"
	case "RhoG":
		l += "\\rho^g"
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
