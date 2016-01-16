// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"os"
	"testing"

	"github.com/cpmech/gofem/mdl/por"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

func Test_msh01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("msh01")

	msh, err := ReadMsh("data", "bh16.msh", 0)
	if err != nil {
		tst.Errorf("test failed:\n%v", err)
		return
	}
	io.Pforan("%v\n", msh)
	io.Pfcyan("lims = [%g, %g, %g, %g, %g, %g]\n", msh.Xmin, msh.Xmax, msh.Ymin, msh.Ymax, msh.Zmin, msh.Zmax)
	chk.Scalar(tst, "xmin", 1e-17, msh.Xmin, 10)
	chk.Scalar(tst, "xmax", 1e-17, msh.Xmax, 14)
	chk.Scalar(tst, "ymin", 1e-17, msh.Ymin, -1)
	chk.Scalar(tst, "ymax", 1e-17, msh.Ymax, 1)

	if chk.Verbose {
		msh.Draw2d(false, true, nil)
		plt.SaveD("/tmp/gofem", "test_msh01.png")
	}
}

func Test_msh02(tst *testing.T) {

	//verbose()
	chk.PrintTitle("msh02")

	msh, err := ReadMsh("data", "ex-beam-joint-comp.msh", 0)
	if err != nil {
		tst.Errorf("test failed:\n%v", err)
		return
	}

	v7 := msh.Verts[7]
	v21 := msh.Verts[21]
	v30 := msh.Verts[30]
	chk.Ints(tst, "vert #  7: SharedBy", v7.SharedBy, []int{0, 1, 2, 3, 4, 5, 6, 7})
	chk.Ints(tst, "vert # 21: SharedBy", v21.SharedBy, []int{4, 5, 6, 7, 8, 9, 10, 11})
	chk.Ints(tst, "vert # 30: SharedBy", v30.SharedBy, []int{8, 9, 10, 11})

	chk.Ints(tst, "vids connected to joint 13", msh.Cells[13].JntConVerts, []int{7, 21})
	chk.Ints(tst, "cids connected to joint 13", msh.Cells[13].JntConCells, []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11})
	chk.Ints(tst, "vids connected to joint 15", msh.Cells[15].JntConVerts, []int{21, 30})
	chk.Ints(tst, "cids connected to joint 15", msh.Cells[15].JntConCells, []int{4, 5, 6, 7, 8, 9, 10, 11})

	io.Pforan("MaxElev = %v\n", msh.MaxElev)
	chk.Scalar(tst, "MaxElev", 1e-17, msh.MaxElev, 1)
}

func Test_msh03(tst *testing.T) {

	//verbose()
	chk.PrintTitle("msh03")

	msh, err := ReadMsh("data", "bjointcomp3d.msh", 0)
	if err != nil {
		tst.Errorf("test failed:\n%v", err)
		return
	}

	io.Pforan(" 0 (hex8)  : IsSolid = %v\n", msh.Cells[0].IsSolid)
	io.Pforan(" 8 (beam)  : IsSolid = %v\n", msh.Cells[8].IsSolid)
	io.Pforan(" 9 (beam)  : IsSolid = %v\n", msh.Cells[9].IsSolid)
	io.Pforan("10 (joint) : IsSolid = %v\n", msh.Cells[10].IsSolid)
	if !msh.Cells[0].IsSolid {
		tst.Errorf("test failed\n")
		return
	}
	if msh.Cells[8].IsSolid {
		tst.Errorf("test failed\n")
		return
	}
	if msh.Cells[9].IsSolid {
		tst.Errorf("test failed\n")
		return
	}
	if msh.Cells[10].IsSolid {
		tst.Errorf("test failed\n")
		return
	}

	io.Pforan("MaxElev = %v\n", msh.MaxElev)
	chk.Scalar(tst, "MaxElev", 1e-17, msh.MaxElev, 1)
}

func Test_mat01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("mat01")

	H := 10.0
	grav := 10.0

	mdb1, err := ReadMat("data", "bh.mat", 2, false, H, grav)
	if err != nil {
		tst.Errorf("cannot read bh.mat\n:%v", err)
		return
	}
	io.Pforan("bh.mat just read:\n%v\n", mdb1)

	fn := "test_bh.mat"
	io.WriteFileSD("/tmp/gofem/inp", fn, mdb1.String())

	mdb2, err := ReadMat("/tmp/gofem/inp/", fn, 2, false, H, grav)
	if err != nil {
		tst.Errorf("cannot read test_bh.mat\n:%v", err)
		return
	}
	io.Pfblue2("\n%v\n", mdb2)
}

func Test_mat02(tst *testing.T) {

	//verbose()
	chk.PrintTitle("mat02")

	H := 10.0
	grav := 10.0

	mdb, err := ReadMat("data", "porous.mat", 2, false, H, grav)
	if err != nil {
		tst.Errorf("cannot read porous.mat\n:%v", err)
		return
	}
	io.Pfblue2("porous.mat just read:\n%v\n", mdb)

	mat := mdb.Get("porous1")
	io.Pforan("Porous = %v\n", mat.Por)
	if chk.Verbose {
		plt.SetForEps(1.2, 400)
		por.PlotLrm(mat.Por, "/tmp/gofem", "fig_mat02_lrm.eps", 20, 101, true, true, true, "", "", "", "", "")
	}
}

func Test_sim01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("sim01")

	sim := ReadSim("data/bh16.sim", "", true, 0)
	if sim == nil {
		tst.Errorf("test failed:\n")
		return
	}
	if chk.Verbose {
		sim.GetInfo(os.Stdout)
		io.Pf("\n")
	}

	io.Pfyel("Ndim    = %v\n", sim.Ndim)
	io.Pfyel("MaxElev = %v\n", sim.MaxElev)

	io.Pfcyan("\nLinSol.Name = %v\n", sim.LinSol.Name)

	chk.IntAssert(sim.Ndim, 2)
	chk.Scalar(tst, "maxElev", 1e-15, sim.MaxElev, 1)
}

func Test_sim02(tst *testing.T) {

	//verbose()
	chk.PrintTitle("sim02")

	sim := ReadSim("data/frees01.sim", "", true, 0)
	if sim == nil {
		tst.Errorf("test failed:\n")
		return
	}
	if chk.Verbose {
		sim.GetInfo(os.Stdout)
		io.Pf("\n")
	}

	io.Pfyel("Ndim    = %v\n", sim.Ndim)
	io.Pfyel("MaxElev = %v\n", sim.MaxElev)

	io.Pf("\n")
	io.Pfyel("liquid: RhoL0 = %v\n", sim.LiqMdl.R0)
	io.Pfyel("liquid: pl0   = %v\n", sim.LiqMdl.P0)
	io.Pfyel("liquid: Cl    = %v\n", sim.LiqMdl.C)
	io.Pfyel("liquid: g     = %v\n", sim.LiqMdl.Grav)
	io.Pfyel("liquid: H     = %v\n", sim.LiqMdl.H)

	io.Pf("\n")
	io.Pf("gas: RhoG0 = %v\n", sim.GasMdl.R0)
	io.Pf("gas: pg0   = %v\n", sim.GasMdl.P0)
	io.Pf("gas: Cg    = %v\n", sim.GasMdl.C)
	io.Pf("gas: g     = %v\n", sim.GasMdl.Grav)
	io.Pf("gas: H     = %v\n", sim.GasMdl.H)

	io.Pf("\n")
	chk.IntAssert(sim.Ndim, 2)
	chk.Scalar(tst, "MaxElev", 1e-15, sim.MaxElev, 10.0)
	io.Pf("\n")
	chk.Scalar(tst, "liq: RhoL0", 1e-15, sim.LiqMdl.R0, 1.0)
	chk.Scalar(tst, "liq: Pl0  ", 1e-15, sim.LiqMdl.P0, 0.0)
	chk.Scalar(tst, "liq: Cl   ", 1e-15, sim.LiqMdl.C, 4.53e-07)
	chk.Scalar(tst, "liq: Grav ", 1e-15, sim.LiqMdl.Grav, 10.0)
	chk.Scalar(tst, "liq: H    ", 1e-15, sim.LiqMdl.H, 10.0)
	io.Pf("\n")
	chk.Scalar(tst, "gas: RhoG0", 1e-15, sim.GasMdl.R0, 0.0012)
	chk.Scalar(tst, "gas: Pg0  ", 1e-15, sim.GasMdl.P0, 0.0)
	chk.Scalar(tst, "gas: Cg   ", 1e-15, sim.GasMdl.C, 1.17e-5)
	chk.Scalar(tst, "gas: Grav ", 1e-15, sim.GasMdl.Grav, 10.0)
	chk.Scalar(tst, "gas: H    ", 1e-15, sim.GasMdl.H, 10.0)

	if chk.Verbose {
		sim.LiqMdl.Plot("/tmp/gofem", "fig_sim02_liq", 21)
		sim.GasMdl.Plot("/tmp/gofem", "fig_sim02_gas", 21)
	}
}
