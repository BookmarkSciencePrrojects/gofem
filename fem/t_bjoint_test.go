// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

func Test_bjoint01a(tst *testing.T) {

	//   Ids and        <10> 21 22 23               LEGEND:
	//   equations:       *                             * beam element
	//                 {4}*                             # beam + joint elements
	//                    *                          [id] id of solids
	//              12 13 * 18 19 20                 {id} id of beams
	//    14 <6>--------<7,9>--------<8> 16          (id) id of joints
	//    15  |           #           |  17          <id> id of nodes
	//        |        {5}#           |             <a,b> id of solid, id of beam node
	//        |    [2]    #    [3]    |
	//        |        (6)#           |
	//        |           #           |
	//     6 <3>--------<4,11>-------<5> 10
	//     7  |         4 | 24        |  11
	//        |         5 | 25        |
	//        |    [0]    | 26 [1]    |
	//        |           |           |
	//        |           |           |
	//     0 <0>---------<1>---------<2> 8
	//     1             2 3             9

	//verbose()
	chk.PrintTitle("bjoint01a. beam joint compatible. static. check")

	// start simulation
	analysis := NewFEM("data/bjointcomp2d01.sim", "", true, false, false, false, chk.Verbose, 0)

	// set stage
	err := analysis.SetStage(0)
	if err != nil {
		tst.Errorf("SetStage failed:\n%v", err)
		return
	}

	// initialise solution vectros
	err = analysis.ZeroStage(0, true)
	if err != nil {
		tst.Errorf("ZeroStage failed:\n%v", err)
		return
	}

	// domain and elements
	dom := analysis.Domains[0]
	solids := []*ElemU{dom.Elems[0].(*ElemU), dom.Elems[1].(*ElemU), dom.Elems[2].(*ElemU), dom.Elems[3].(*ElemU)}
	beams := []*Beam{dom.Elems[4].(*Beam), dom.Elems[5].(*Beam)}
	joint := dom.Elems[6].(*BjointComp)

	// nodes and elements
	chk.IntAssert(len(dom.Nodes), 12)
	chk.IntAssert(len(dom.Elems), 7)

	// check dofs
	nods_with_rz := map[int]bool{11: true, 9: true, 10: true}
	for _, nod := range dom.Nodes {
		if nods_with_rz[nod.Vert.Id] {
			chk.IntAssert(len(nod.Dofs), 3)
			chk.StrAssert(nod.Dofs[0].Key, "ux")
			chk.StrAssert(nod.Dofs[1].Key, "uy")
			chk.StrAssert(nod.Dofs[2].Key, "rz")
		} else {
			chk.IntAssert(len(nod.Dofs), 2)
			chk.StrAssert(nod.Dofs[0].Key, "ux")
			chk.StrAssert(nod.Dofs[1].Key, "uy")
		}
	}

	// check equations
	nids, eqs := get_nids_eqs(dom)
	chk.Ints(tst, "eqs", eqs, utl.IntRange(27))
	chk.Ints(tst, "nids", nids, []int{0, 1, 4, 3, 2, 5, 7, 6, 8, 9, 10, 11})

	// check Umaps
	chk.Ints(tst, "solid 0 Umap", solids[0].Umap, []int{0, 1, 2, 3, 4, 5, 6, 7})
	chk.Ints(tst, "solid 1 Umap", solids[1].Umap, []int{2, 3, 8, 9, 10, 11, 4, 5})
	chk.Ints(tst, "solid 2 Umap", solids[2].Umap, []int{6, 7, 4, 5, 12, 13, 14, 15})
	chk.Ints(tst, "solid 3 Umap", solids[3].Umap, []int{4, 5, 10, 11, 16, 17, 12, 13})
	chk.Ints(tst, "beam  0 Umap", beams[0].Umap, []int{18, 19, 20, 21, 22, 23})
	chk.Ints(tst, "beam  1 Umap", beams[1].Umap, []int{24, 25, 26, 18, 19, 20})
	chk.Ints(tst, "joint SldLocVid", joint.SldLocVid, []int{1, 2})
	chk.Ints(tst, "joint 0 LinUmap", joint.LinUmap, []int{24, 25, 18, 19})
	chk.Ints(tst, "joint 0 SldUmap", joint.SldUmap, []int{4, 5, 12, 13})

	// intial values @ integration points
	io.Pforan("initial values @ integration points\n")
	for _, ele := range solids {
		for idx, ip := range ele.IpsElem {
			s := ele.States[idx]
			y := ele.Cell.Shp.IpRealCoords(ele.X, ip)[1]
			chk.AnaNum(tst, io.Sf("σ(y=%11.8f)", y), 1e-17, s.Sig[0], -1, chk.Verbose)
		}
	}

	// joint
	io.Pforan("e0 = %v\n", joint.Lin.e0)
	io.Pforan("e1 = %v\n", joint.Lin.e1)
	io.Pforan("e2 = %v\n", joint.Lin.e2)

	// pressure
	τ := joint.States[0].Sig
	q1 := joint.States[0].Phi[0]
	q2 := joint.States[0].Phi[1]
	io.Pfyel("τ=%g q1=%g q2=%g\n", τ, q1, q2)
	chk.Scalar(tst, "τ", 1e-17, τ, 0)
	chk.Scalar(tst, "q1", 1e-15, q1, 0)
	chk.Scalar(tst, "q2", 1e-15, q2, 0)

	// contact force vector
	joint.AddToRhs(dom.Fb, dom.Sol)
	io.Pforan("fC = %v\n", joint.fC)
	io.Pfblue2("fb = %v\n", dom.Fb)

	ele2 := dom.Elems[2].(*ElemU)
	io.Pfgreen("sld: sig = %v\n", ele2.States[0].Sig)
}

func Test_bjoint01b(tst *testing.T) {

	//verbose()
	chk.PrintTitle("bjoint01b. beam joint compatible. static. run")

	// start simulation
	analysis := NewFEM("data/bjointcomp2d01.sim", "", true, true, false, false, chk.Verbose, 0)

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
	}
}

func Test_bjoint02(tst *testing.T) {

	//verbose()
	chk.PrintTitle("bjoint02. beam joint compatible. pull-out. run")

	// start simulation
	analysis := NewFEM("data/bjointcomp2d02.sim", "", true, true, false, false, chk.Verbose, 0)

	// for debugging Kb
	if true {
		bjointcomp_DebugKb(analysis, &testKb{tst: tst, verb: chk.Verbose, eid: 6,
			ni: -1, nj: -1, itmin: 1, itmax: -1, tmin: -1, tmax: -1, tol: 1e-10,
		})
	}

	// run simulation
	err := analysis.Run()
	if err != nil {
		tst.Errorf("Run failed:\n%v", err)
	}

	// displacements @ tip
	dom := analysis.Domains[0]
	tip := dom.Msh.VertTag2verts[-102][0].Id
	nod := dom.Nodes[tip]
	eq := nod.GetEq("uy")
	io.Pforan("uy @ tip (%d,%d) = %v\n", -102, tip, dom.Sol.Y[eq])

	// stresses along embedded beam
	emb := dom.Elems[6].(*BjointComp)
	io.Pfblue2("τ = %v\n", emb.States[0].Sig)
}
