// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import "github.com/cpmech/gosl/io"

func (o Point) String() string {
	l := io.Sf("\"Vid:%d\", \"IpId:%d\", \"dist\":%g, \"X\":[%g,%g", o.Vid, o.IpId, o.Dist, o.X[0], o.X[1])
	if len(o.X) == 3 {
		l += io.Sf(",%g", o.X[2])
	}
	l += "], \"vals\":["
	first := true
	for key, vals := range o.Vals {
		if !first {
			l += ", "
		}
		l += io.Sf("{\"key\":%q, nVals=%d}", key, len(vals))
		first = false
	}
	l += "]}"
	return l
}

func (o Points) String() string {
	l := "[\n"
	for i, q := range o {
		if i > 0 {
			l += ",\n"
		}
		l += io.Sf("    %v", q)
	}
	l += "\n  ]"
	return l
}

func (o ResultsMap) String() string {
	l := "{\n"
	first := true
	for key, val := range o {
		if !first {
			l += ",\n"
		}
		l += io.Sf("  %q : %v", key, val)
		first = false
	}
	if len(o) > 0 {
		l += "\n"
	}
	l += "}"
	return l
}
