// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/io"
)

func main() {

	// filename
	filename, _ := io.ArgToFilename(0, "patch", ".sim", true)

	// results
	out.Start(filename, 0, 0)
	out.Define("p0", out.N{-100}) // -100=tag
	out.Define("p1", out.N{-101}) // -101=tag
	out.Define("p2", out.N{-102}) // -102=tag
	out.Define("p3", out.N{-103}) // -103=tag
	out.Define("p4", out.N{-104}) // -104=tag
	out.Define("p5", out.N{-105}) // -105=tag
	out.Define("p6", out.N{-106}) // -106=tag
	out.Define("p7", out.N{-107}) // -107=tag
	out.LoadResults(nil)
	
	time := 1
	ux0 := out.GetRes("ux", "p0", -1)[time]
	ux1 := out.GetRes("ux", "p1", -1)[time]
	ux2 := out.GetRes("ux", "p2", -1)[time]
	ux3 := out.GetRes("ux", "p3", -1)[time]
	ux4 := out.GetRes("ux", "p4", -1)[time]
	ux5 := out.GetRes("ux", "p5", -1)[time]
	ux6 := out.GetRes("ux", "p6", -1)[time]
	ux7 := out.GetRes("ux", "p7", -1)[time]
	uy0 := out.GetRes("uy", "p0", -1)[time]
	uy1 := out.GetRes("uy", "p1", -1)[time]
	uy2 := out.GetRes("uy", "p2", -1)[time]
	uy3 := out.GetRes("uy", "p3", -1)[time]
	uy4 := out.GetRes("uy", "p4", -1)[time]
	uy5 := out.GetRes("uy", "p5", -1)[time]
	uy6 := out.GetRes("uy", "p6", -1)[time]
	uy7 := out.GetRes("uy", "p7", -1)[time]

	// show FEM results
	io.Pforan("ux_p0 = %v\n", ux0)
	io.Pforan("ux_p1 = %v\n", ux1)
	io.Pforan("ux_p2 = %v\n", ux2)
	io.Pforan("ux_p3 = %v\n", ux3)
	io.Pforan("ux_p4 = %v\n", ux4)
	io.Pforan("ux_p5 = %v\n", ux5)
	io.Pforan("ux_p6 = %v\n", ux6)
	io.Pforan("ux_p7 = %v\n", ux7)
	io.Pforan("uy_p0 = %v\n", uy0)
	io.Pforan("uy_p1 = %v\n", uy1)
	io.Pforan("uy_p2 = %v\n", uy2)
	io.Pforan("uy_p3 = %v\n", uy3)
	io.Pforan("uy_p4 = %v\n", uy4)
	io.Pforan("uy_p5 = %v\n", uy5)
	io.Pforan("uy_p6 = %v\n", uy6)
	io.Pforan("uy_p7 = %v\n", uy7)
}
