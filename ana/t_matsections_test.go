// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package ana

import (
	"math"
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func Test_sections01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("sections01. typical cross-sections")

	var rect CrossSection
	b, h := 4.0, 6.0
	rect.Init("rectangle", "in", b, h, 0, 0, 0)
	io.Pforan("4 x 6 rectangle:\n%v\n", rect.GetMatString("%g"))
	chk.Scalar(tst, "rect: A  ", 1e-17, rect.A, 24.0)
	chk.Scalar(tst, "rect: I22", 1e-17, rect.I22, 72.0)
	chk.Scalar(tst, "rect: I11", 1e-17, rect.I11, 32.0)
	chk.Scalar(tst, "rect: Jtt", 1e-11, rect.Jtt, 75.1249382716)

	b, h = 4.0, 4.0
	rect.Init("rectangle", "in", b, h, 0, 0, 0)
	io.Pforan("\n4 x 4 rectangle:\n%v\n", rect.GetMatString("%g"))
	chk.Scalar(tst, "rect: A  ", 1e-17, rect.A, 16.0)
	chk.Scalar(tst, "rect: I22", 1e-13, rect.I22, 21.3333333333333)
	chk.Scalar(tst, "rect: I11", 1e-13, rect.I11, 21.3333333333333)
	chk.Scalar(tst, "rect: Jtt", 1e-17, rect.Jtt, 36.0)

	var ibeam CrossSection
	b, h = 4.0, 6.0
	tf, tw := 0.5, 0.3
	ibeam.Init("I-beam", "in", b, h, tf, tw, 0)
	io.Pforan("\n4 x 6 I-beam:\n%v\n", ibeam.GetMatString("%g"))
	chk.Scalar(tst, "I-beam: A  ", 1e-17, ibeam.A, 5.5)
	chk.Scalar(tst, "I-beam: I22", 1e-10, ibeam.I22, 33.4583333333)
	chk.Scalar(tst, "I-beam: I11", 1e-10, ibeam.I11, 5.3445833333)
	chk.Scalar(tst, "I-beam: Jtt", 1e-10, ibeam.Jtt, 0.3783333333)

	var circle CrossSection
	r := 1.0
	circle.Init("circle", "m", 0, 0, 0, 0, r)
	io.Pforan("\nr=1 circle:\n%v\n", circle.GetMatString("%g"))
	chk.Scalar(tst, "circle: A  ", 1e-17, circle.A, math.Pi)
	chk.Scalar(tst, "circle: I22", 1e-10, circle.I22, 0.7853981634)
	chk.Scalar(tst, "circle: I11", 1e-10, circle.I11, 0.7853981634)
	chk.Scalar(tst, "circle: Jtt", 1e-11, circle.Jtt, 1.5707963268)
}

func Test_materials01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("materials01. reference materials parameters")

	var rect CrossSection
	b, h := 0.2, 0.3
	rect.Init("rectangle", "m", b, h, 0, 0, 0)

	var mat Material
	mat.Init("steel", "MPa")
	io.Pforan("%v\n", mat.GetMatString("steel", "", "%e", &rect))
}
