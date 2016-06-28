// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"encoding/json"
	"math"
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/num"
)

// T_iteration testing: iteration results
type T_iteration struct {
	It     int     // iteration number
	ResRel float64 // relative residual
	Resid  float64 // absolute residual
}

// T_results testing: results
type T_results struct {
	Status     string        // status message
	LoadFactor float64       // load factor
	Iterations []T_iteration // iterations data
	Kmats      [][][]float64 // [nele][nu][nu] all stiffness matrices
	Disp       [][]float64   // [nnod][ndim] displacements at nodes
	DispMult   float64       // displacements multiplier
	Note       string        // note about number of integration points
	Stresses   [][][]float64 // [nele][nip][nsig] all stresses @ all ips
	// Stresses: examples:
	//   2D solid: [sx, sy, sxy, sz]
	//   3D beam:  [M22, M11, Mtt, Vs, Vr]
}

// T_results_set is a set of comparison results
type T_results_set []*T_results

// testing_compare_results_u compares results with u-formulation
func TestingCompareResultsU(tst *testing.T, simfilepath, cmpfname, alias string, tolK, tolu, tols float64, skipK, verbose bool, extraAfterSetStage func(dom *Domain)) {

	// flag
	do_check_stresses := true

	// FEM structure
	fem := NewFEM(simfilepath, alias, false, false, true, false, verbose, 0)

	// set stage
	err := fem.SetStage(0)
	if err != nil {
		chk.Panic("cannot set stage:\n%v", err)
	}

	// extra settings
	if extraAfterSetStage != nil {
		extraAfterSetStage(fem.Domains[0])
	}

	// zero solution
	err = fem.ZeroStage(0, true)
	if err != nil {
		chk.Panic("cannot zero stage data:\n%v", err)
	}

	// read file with comparison results
	buf, err := io.ReadFile(cmpfname)
	if err != nil {
		tst.Errorf("TestingCompareResultsU: ReadFile failed\n")
		return
	}

	// unmarshal json
	var cmp_set T_results_set
	err = json.Unmarshal(buf, &cmp_set)
	if err != nil {
		tst.Errorf("TestingCompareResultsU: Unmarshal failed\n")
		return
	}

	// run comparisons
	dom := fem.Domains[0]
	dmult := 1.0
	for idx, cmp := range cmp_set {

		// displacements multiplier
		if idx == 0 && math.Abs(cmp.DispMult) > 1e-10 {
			dmult = cmp.DispMult
		}

		// time index
		tidx := idx + 1
		if verbose {
			io.PfYel("\n\ntidx = %d . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .\n", tidx)
		}

		// load gofem results
		err = dom.Read(fem.Summary, tidx, 0, true)
		if err != nil {
			chk.Panic("cannot read 'gofem' results:\n%v", err)
		}
		if verbose {
			io.Pfyel("time = %v\n", dom.Sol.T)
		}

		// check K matrices
		if !skipK {
			if verbose {
				io.Pfgreen(". . . checking K matrices . . .\n")
			}
			for eid, Ksg := range cmp.Kmats {
				if e, ok := dom.Elems[eid].(*ElemU); ok {
					err = e.AddToKb(dom.Kb, dom.Sol, true)
					if err != nil {
						chk.Panic("TestingCompareResultsU: AddToKb failed\n")
					}
					chk.Matrix(tst, io.Sf("K%d", eid), tolK, e.K, Ksg)
				}
				if e, ok := dom.Elems[eid].(*Beam); ok {
					err = e.AddToKb(dom.Kb, dom.Sol, true)
					if err != nil {
						chk.Panic("TestingCompareResultsU: AddToKb failed\n")
					}
					chk.Matrix(tst, io.Sf("K%d", eid), tolK, e.K, Ksg)
				}
				if e, ok := dom.Elems[eid].(*Rod); ok {
					err = e.AddToKb(dom.Kb, dom.Sol, true)
					if err != nil {
						chk.Panic("TestingCompareResultsU: AddToKb failed\n")
					}
					chk.Matrix(tst, io.Sf("K%d", eid), tolK, e.K, Ksg)
				}
				if e, ok := dom.Elems[eid].(*ElastRod); ok {
					err = e.AddToKb(dom.Kb, dom.Sol, true)
					if err != nil {
						chk.Panic("TestingCompareResultsU: AddToKb failed\n")
					}
					chk.Matrix(tst, io.Sf("K%d", eid), tolK, e.K, Ksg)
				}
			}
		}

		// check displacements
		if verbose {
			io.Pfgreen(". . . checking displacements . . .\n")
		}
		for nid, usg := range cmp.Disp {
			ix := dom.Vid2node[nid].Dofs[0].Eq
			iy := dom.Vid2node[nid].Dofs[1].Eq
			chk.AnaNum(tst, "ux", tolu, dom.Sol.Y[ix], usg[0]*dmult, verbose)
			chk.AnaNum(tst, "uy", tolu, dom.Sol.Y[iy], usg[1]*dmult, verbose)
			if len(usg) > 2 {
				iz := dom.Vid2node[nid].Dofs[2].Eq
				chk.AnaNum(tst, "uz", tolu, dom.Sol.Y[iz], usg[2]*dmult, verbose)
				for j := 3; j < len(usg); j++ {
					idx := dom.Vid2node[nid].Dofs[j].Eq
					chk.AnaNum(tst, io.Sf("u%d", j), tolu, dom.Sol.Y[idx], usg[j]*dmult, verbose)
				}
			}
		}

		// check stresses
		if do_check_stresses {
			if verbose {
				io.Pfgreen(". . . checking stresses (or moments/shear forces). . .\n")
			}
			for eid, sig := range cmp.Stresses {
				if verbose {
					io.Pforan("eid = %d\n", eid)
				}
				if e, ok := dom.Cid2elem[eid].(*ElemU); ok {
					for ip, val := range sig {
						if verbose {
							io.Pfgrey2("ip = %d @ %v\n", ip, e.Cell.Shp.IpRealCoords(e.X, e.IpsElem[ip]))
						}
						σ := e.States[ip].Sig
						if len(val) == 6 {
							chk.AnaNum(tst, "sx ", tols, σ[0], val[0], verbose)
							chk.AnaNum(tst, "sy ", tols, σ[1], val[1], verbose)
						} else {
							chk.AnaNum(tst, "sx ", tols, σ[0], val[0], verbose)
							chk.AnaNum(tst, "sy ", tols, σ[1], val[1], verbose)
							chk.AnaNum(tst, "sxy", tols, σ[3]/SQ2, val[2], verbose)
							if len(val) > 3 { // sx, sy, sxy, sz
								chk.AnaNum(tst, "sz ", tols, σ[2], val[3], verbose)
							}
						}
					}
				}
				if e, ok := dom.Cid2elem[eid].(*Rod); ok {
					for ip, val := range sig {
						if verbose {
							io.Pfgrey2("ip = %d\n", ip)
						}
						σ := e.States[ip].Sig
						chk.AnaNum(tst, "sig", tols, σ, val[0], verbose)
					}
				}
				if e, ok := dom.Cid2elem[eid].(*ElastRod); ok {
					res := NewIpsMap()
					e.OutIpVals(res, dom.Sol)
					for ip, val := range sig {
						if verbose {
							io.Pfgrey2("ip = %d\n", ip)
						}
						σ := res.Get("sig", 0)
						chk.AnaNum(tst, "sig", tols, σ, val[0], verbose)
					}
				}
				if e, ok := dom.Cid2elem[eid].(*Beam); ok {
					res := NewIpsMap()
					e.OutIpVals(res, dom.Sol)
					for st, val := range sig {
						if verbose {
							io.Pfgrey2("station = %d\n", st)
						}
						if e.Ndim == 2 {
							M22 := res.Get("M22", st)
							chk.AnaNum(tst, "M22", tols, M22, val[0], verbose)
						} else {
							M22, M11, Mtt := res.Get("M22", st), res.Get("M11", st), res.Get("T00", st)
							chk.AnaNum(tst, "M22", tols, M22, val[0], verbose)
							chk.AnaNum(tst, "M11", tols, M11, val[1], verbose)
							chk.AnaNum(tst, "Mtt", tols, Mtt, val[2], verbose)
						}
					}
				}
			}
		}
	}
}

// testKb helps on checking Kb matrices
type testKb struct {

	// input (must)
	tst          *testing.T // testing structure
	eid          int        // element id
	tol          float64    // tolerance to compare K's
	tol2         float64    // another tolerance to compare K's
	step         float64    // step for finite differences method
	verb         bool       // verbose: show results
	ni, nj       int        // number of i and j components of K to be tested; -1 means all K components
	itmin, itmax int        // limits to consider test; -1 means all iterations
	tmin, tmax   float64    // limits to consider test; -1 means all times

	// derived
	it    int       // current iteration
	t     float64   // current time
	Fbtmp []float64 // auxiliary array
	Ybkp  []float64 // auxiliary array
	ΔYbkp []float64 // auxiliary array
	Yold  []float64 // auxiliary array
}

// p_DebugKb defines a global function to debug Kb for p-elements
func p_DebugKb(fem *FEM, o *testKb) {
	fem.DebugKb = func(d *Domain, it int) {

		elem := d.Elems[o.eid]
		if e, ok := elem.(*ElemP); ok {

			// skip?
			o.it = it
			o.t = d.Sol.T
			if o.skip() {
				return
			}

			// backup and restore upon exit
			o.aux_backup(d)
			defer func() { o.aux_restore(d) }()

			// check
			o.check("Kpp", d, e, e.Pmap, e.Pmap, e.Kpp, o.tol)
			o.check("Kpf", d, e, e.Pmap, e.Fmap, e.Kpf, o.tol)
			o.check("Kfp", d, e, e.Fmap, e.Pmap, e.Kfp, o.tol)
			o.check("Kff", d, e, e.Fmap, e.Fmap, e.Kff, o.tol)
		}
	}
}

// pp_DebugKb defines a global function to debug Kb for pp-elements
func pp_DebugKb(fem *FEM, o *testKb) {
	fem.DebugKb = func(d *Domain, it int) {

		elem := d.Elems[o.eid]
		if e, ok := elem.(*ElemPP); ok {

			// skip?
			o.it = it
			o.t = d.Sol.T
			if o.skip() {
				return
			}

			// backup and restore upon exit
			o.aux_backup(d)
			defer func() { o.aux_restore(d) }()

			// check
			o.check("Kll", d, e, e.Plmap, e.Plmap, e.Kll, o.tol)
			o.check("Klg", d, e, e.Plmap, e.Pgmap, e.Klg, o.tol)
			o.check("Kgl", d, e, e.Pgmap, e.Plmap, e.Kgl, o.tol)
			o.check("Kgg", d, e, e.Pgmap, e.Pgmap, e.Kgg, o.tol)
			o.check("Klf", d, e, e.Plmap, e.Flmap, e.Klf, o.tol)
			o.check("Kfl", d, e, e.Flmap, e.Plmap, e.Kfl, o.tol)
			o.check("Kff", d, e, e.Flmap, e.Flmap, e.Kff, o.tol)
		}
	}
}

// u_DebugKb defines a global function to debug Kb for u-elements
func u_DebugKb(fem *FEM, o *testKb) {
	fem.DebugKb = func(d *Domain, it int) {

		elem := d.Elems[o.eid]
		if e, ok := elem.(*ElemU); ok {

			// skip?
			o.it = it
			o.t = d.Sol.T
			if o.skip() {
				return
			}

			// backup and restore upon exit
			o.aux_backup(d)
			defer func() { o.aux_restore(d) }()

			// check
			if e.HasContact {
				o.check("Kqq", d, e, e.Qmap, e.Qmap, e.Kqq, o.tol)
				o.check("Kqu", d, e, e.Qmap, e.Umap, e.Kqu, o.tol)
				o.check("Kuq", d, e, e.Umap, e.Qmap, e.Kuq, o.tol)
			}
			o.check("K", d, e, e.Umap, e.Umap, e.K, o.tol)
		}
	}
	return
}

// up_DebugKb defines a global function to debug Kb for up-elements
func up_DebugKb(fem *FEM, o *testKb) {
	fem.DebugKb = func(d *Domain, it int) {

		elem := d.Elems[o.eid]
		if e, ok := elem.(*ElemUP); ok {

			// skip?
			o.it = it
			o.t = d.Sol.T
			if o.skip() {
				return
			}

			// backup and restore upon exit
			o.aux_backup(d)
			defer func() { o.aux_restore(d) }()

			// check
			o.check("Kuu", d, e, e.U.Umap, e.U.Umap, e.U.K, o.tol)
			o.check("Kup", d, e, e.U.Umap, e.P.Pmap, e.Kup, o.tol2)
			o.check("Kpu", d, e, e.P.Pmap, e.U.Umap, e.Kpu, o.tol)
			o.check("Kpp", d, e, e.P.Pmap, e.P.Pmap, e.P.Kpp, o.tol)
			o.check("Kpf", d, e, e.P.Pmap, e.P.Fmap, e.P.Kpf, o.tol)
			o.check("Kfp", d, e, e.P.Fmap, e.P.Pmap, e.P.Kfp, o.tol)
			o.check("Kff", d, e, e.P.Fmap, e.P.Fmap, e.P.Kff, o.tol)
		}
	}
	return
}

// upp_DebugKb defines a global function to debug Kb for upp-elements
func upp_DebugKb(fem *FEM, o *testKb) {
	fem.DebugKb = func(d *Domain, it int) {

		elem := d.Elems[o.eid]
		if e, ok := elem.(*ElemUPP); ok {

			// skip?
			o.it = it
			o.t = d.Sol.T
			if o.skip() {
				return
			}

			// backup and restore upon exit
			o.aux_backup(d)
			defer func() { o.aux_restore(d) }()

			// {pl,pg,fl} versus {pl,pg,fl}
			o.check("Kll", d, e, e.P.Plmap, e.P.Plmap, e.P.Kll, o.tol)
			o.check("Klg", d, e, e.P.Plmap, e.P.Pgmap, e.P.Klg, o.tol)
			o.check("Kgl", d, e, e.P.Pgmap, e.P.Plmap, e.P.Kgl, o.tol)
			o.check("Kgg", d, e, e.P.Pgmap, e.P.Pgmap, e.P.Kgg, o.tol)
			o.check("Klf", d, e, e.P.Plmap, e.P.Flmap, e.P.Klf, o.tol)
			o.check("Kfl", d, e, e.P.Flmap, e.P.Plmap, e.P.Kfl, o.tol)
			o.check("Kff", d, e, e.P.Flmap, e.P.Flmap, e.P.Kff, o.tol)

			// {pl,pg,u} versus {pl,pg,u}
			o.check("Kul", d, e, e.U.Umap, e.P.Plmap, e.Kul, o.tol2)
			o.check("Kug", d, e, e.U.Umap, e.P.Pgmap, e.Kug, o.tol2)
			o.check("Klu", d, e, e.P.Plmap, e.U.Umap, e.Klu, o.tol)
			o.check("Kgu", d, e, e.P.Pgmap, e.U.Umap, e.Kgu, o.tol)

			// {u} versus {u}
			o.check("Kuu", d, e, e.U.Umap, e.U.Umap, e.U.K, o.tol)
		}
	}
	return
}

// rjoint_DebugKb defines a global function to debug Kb for rjoint-elements
func rjoint_DebugKb(fem *FEM, o *testKb) {
	fem.DebugKb = func(d *Domain, it int) {

		elem := d.Elems[o.eid]
		if e, ok := elem.(*Rjoint); ok {

			// skip?
			o.it = it
			o.t = d.Sol.T
			if o.skip() {
				return
			}

			// backup and restore upon exit
			o.aux_backup(d)
			defer func() { o.aux_restore(d) }()

			// check
			o.check("Krr", d, e, e.Rod.Umap, e.Rod.Umap, e.Krr, o.tol)
			o.check("Krs", d, e, e.Rod.Umap, e.Sld.Umap, e.Krs, o.tol)
			o.check("Ksr", d, e, e.Sld.Umap, e.Rod.Umap, e.Ksr, o.tol)
			o.check("Kss", d, e, e.Sld.Umap, e.Sld.Umap, e.Kss, o.tol)
		} else {
			io.Pfred("warning: eid=%d does not correspond to Rjoint element\n", o.eid)
		}
	}
	return
}

// bjointcomp_DebugKb defines a global function to debug Kb for bjoint-compatible elements
func bjointcomp_DebugKb(fem *FEM, o *testKb) {
	fem.DebugKb = func(d *Domain, it int) {

		elem := d.Elems[o.eid]
		if e, ok := elem.(*BjointComp); ok {

			// skip?
			o.it = it
			o.t = d.Sol.T
			if o.skip() {
				return
			}

			// backup and restore upon exit
			o.aux_backup(d)
			defer func() { o.aux_restore(d) }()

			// check
			o.check("Kll", d, e, e.LinUmap, e.LinUmap, e.Kll, o.tol)
			o.check("Kls", d, e, e.LinUmap, e.SldUmap, e.Kls, o.tol)
			o.check("Ksl", d, e, e.SldUmap, e.LinUmap, e.Ksl, o.tol)
			o.check("Kss", d, e, e.SldUmap, e.SldUmap, e.Kss, o.tol)
		} else {
			io.Pfred("warning: eid=%d does not correspond to BjointComp element\n", o.eid)
		}
	}
	return
}

// skip skips test based on it and/or t
func (o testKb) skip() bool {
	if o.itmin >= 0 {
		if o.it < o.itmin {
			return true // skip
		}
	}
	if o.itmax >= 0 {
		if o.it > o.itmax {
			return true // skip
		}
	}
	if o.tmin >= 0 {
		if o.t < o.tmin {
			return true // skip
		}
	}
	if o.tmax >= 0 {
		if o.t > o.tmax {
			return true // skip
		}
	}
	if o.verb {
		io.PfYel("\nit=%2d t=%v\n", o.it, o.t)
	}
	return false
}

// aux_arrays generates auxiliary arrays
func (o *testKb) aux_backup(d *Domain) {
	o.Fbtmp = make([]float64, d.Ny)
	o.Yold = make([]float64, d.Ny)
	o.Ybkp = make([]float64, d.Ny)
	o.ΔYbkp = make([]float64, d.Ny)
	for i := 0; i < d.Ny; i++ {
		o.Yold[i] = d.Sol.Y[i] - d.Sol.ΔY[i]
		o.Ybkp[i] = d.Sol.Y[i]
		o.ΔYbkp[i] = d.Sol.ΔY[i]
	}
	for _, eivs := range d.ElemIntvars {
		eivs.BackupIvs(true) // true => aux
	}
}

func (o *testKb) aux_restore(d *Domain) {
	for i := 0; i < d.Ny; i++ {
		d.Sol.Y[i] = o.Ybkp[i]
		d.Sol.ΔY[i] = o.ΔYbkp[i]
	}
	for _, eivs := range d.ElemIntvars {
		eivs.RestoreIvs(true) // true => aux
	}
}

// check performs the checking of Kb using numerical derivatives
func (o *testKb) check(label string, d *Domain, e Elem, Imap, Jmap []int, Kana [][]float64, tol float64) {
	var imap, jmap []int
	if o.ni < 0 {
		imap = Imap
	} else {
		if o.ni <= len(Imap) {
			imap = Imap[:o.ni]
		}
	}
	if o.nj < 0 {
		jmap = Jmap
	} else {
		if o.nj <= len(Jmap) {
			jmap = Jmap[:o.nj]
		}
	}
	if o.step < 1e-14 {
		o.step = 1e-6
	}
	//derivfcn := num.DerivForward
	//derivfcn := num.DerivBackward
	derivfcn := num.DerivCentral
	var tmp float64
	for i, I := range imap {
		for j, J := range jmap {
			dnum, _ := derivfcn(func(x float64, args ...interface{}) (res float64) {
				tmp, d.Sol.Y[J] = d.Sol.Y[J], x
				for k := 0; k < d.Ny; k++ {
					o.Fbtmp[k] = 0
					d.Sol.ΔY[k] = d.Sol.Y[k] - o.Yold[k]
				}
				for _, eivs := range d.ElemIntvars {
					eivs.RestoreIvs(false)
				}
				err := d.UpdateElems()
				if err != nil {
					chk.Panic("testing: check: cannot update elements")
				}
				e.AddToRhs(o.Fbtmp, d.Sol)
				res = -o.Fbtmp[I]
				d.Sol.Y[J] = tmp
				return res
			}, d.Sol.Y[J], o.step)
			chk.AnaNum(o.tst, io.Sf(label+"%3d%3d", i, j), tol, Kana[i][j], dnum, o.verb)
		}
	}
}
