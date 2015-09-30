// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"encoding/json"
	"math"
	"testing"

	"github.com/cpmech/gofem/mporous"
	"github.com/cpmech/gofem/msolid"
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
func TestingCompareResultsU(tst *testing.T, simfilepath, cmpfname, alias string, tolK, tolu, tols float64, skipK, verbose bool) {

	// flag
	do_check_stresses := true

	// FEM structure
	fem := NewFEM(simfilepath, alias, false, false, true, false, verbose, 0)

	// set stage
	err := fem.SetStage(0)
	if err != nil {
		chk.Panic("cannot set stage:\n%v", err)
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
							io.Pfgrey2("ip = %d\n", ip)
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
					dat := e.OutIpsData()
					for ip, val := range sig {
						if verbose {
							io.Pfgrey2("ip = %d\n", ip)
						}
						res := dat[ip].Calc(dom.Sol)
						σ := res["sig"]
						chk.AnaNum(tst, "sig", tols, σ, val[0], verbose)
					}
				}
				if e, ok := dom.Cid2elem[eid].(*Beam); ok {
					dat := e.OutIpsData()
					if len(dat) != len(sig) {
						tst.Errorf("number of stations in cmp file is different than the number of stations in element. %d != %d", len(dat), len(sig))
						return
					}
					for station, val := range sig {
						if verbose {
							io.Pfgrey2("station = %d\n", station)
						}
						res := dat[station].Calc(dom.Sol)
						if e.Ndim == 3 {
							M22, M11, Mtt := res["M22"], res["M11"], res["Mtt"]
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
	verb         bool       // verbose: show results
	ni, nj       int        // number of i and j components of K to be tested; -1 means all K components
	itmin, itmax int        // limits to consider test; -1 means all iterations
	tmin, tmax   float64    // limits to consider test; -1 means all times

	// derived
	it    int       // current iteration
	t     float64   // current time
	Fbtmp []float64 // auxiliary array
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

			// copy states and solution
			nip := len(e.IpsElem)
			states := make([]*mporous.State, nip)
			statesBkp := make([]*mporous.State, nip)
			for i := 0; i < nip; i++ {
				states[i] = e.States[i].GetCopy()
				statesBkp[i] = e.StatesBkp[i].GetCopy()
			}
			o.aux_arrays(d)

			// make sure to restore states and solution
			defer func() {
				for i := 0; i < nip; i++ {
					e.States[i].Set(states[i])
					e.StatesBkp[i].Set(statesBkp[i])
				}
				copy(d.Sol.ΔY, o.ΔYbkp)
			}()

			// define restore function
			restore := func() {
				if it == 0 {
					for k := 0; k < nip; k++ {
						e.States[k].Set(states[k])
					}
					return
				}
				for k := 0; k < nip; k++ {
					e.States[k].Set(statesBkp[k])
				}
			}

			// check
			o.check("Kpp", d, e, e.Pmap, e.Pmap, e.Kpp, restore)
			o.check("Kpf", d, e, e.Pmap, e.Fmap, e.Kpf, restore)
			o.check("Kfp", d, e, e.Fmap, e.Pmap, e.Kfp, restore)
			o.check("Kff", d, e, e.Fmap, e.Fmap, e.Kff, restore)
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

			// copy states and solution
			nip := len(e.IpsElem)
			states := make([]*msolid.State, nip)
			statesBkp := make([]*msolid.State, nip)
			for i := 0; i < nip; i++ {
				states[i] = e.States[i].GetCopy()
				statesBkp[i] = e.StatesBkp[i].GetCopy()
			}
			o.aux_arrays(d)

			// make sure to restore states and solution
			defer func() {
				for i := 0; i < nip; i++ {
					e.States[i].Set(states[i])
					e.StatesBkp[i].Set(statesBkp[i])
				}
				copy(d.Sol.ΔY, o.ΔYbkp)
			}()

			// define restore function
			restore := func() {
				if it == 0 {
					for k := 0; k < nip; k++ {
						e.States[k].Set(states[k])
					}
					return
				}
				for k := 0; k < nip; k++ {
					e.States[k].Set(statesBkp[k])
				}
			}

			// check
			if e.HasContact {
				o.check("Kqq", d, e, e.Qmap, e.Qmap, e.Kqq, restore)
				o.check("Kqu", d, e, e.Qmap, e.Umap, e.Kqu, restore)
				o.check("Kuq", d, e, e.Umap, e.Qmap, e.Kuq, restore)
			}
			o.check("K", d, e, e.Umap, e.Umap, e.K, restore)
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

			// copy states and solution
			nip := len(e.U.IpsElem)
			u_states := make([]*msolid.State, nip)
			p_states := make([]*mporous.State, nip)
			u_statesBkp := make([]*msolid.State, nip)
			p_statesBkp := make([]*mporous.State, nip)
			for i := 0; i < nip; i++ {
				u_states[i] = e.U.States[i].GetCopy()
				p_states[i] = e.P.States[i].GetCopy()
				u_statesBkp[i] = e.U.StatesBkp[i].GetCopy()
				p_statesBkp[i] = e.P.StatesBkp[i].GetCopy()
			}
			o.aux_arrays(d)

			// make sure to restore states and solution
			defer func() {
				for i := 0; i < nip; i++ {
					e.U.States[i].Set(u_states[i])
					e.P.States[i].Set(p_states[i])
					e.U.StatesBkp[i].Set(u_statesBkp[i])
					e.P.StatesBkp[i].Set(p_statesBkp[i])
				}
				copy(d.Sol.ΔY, o.ΔYbkp)
			}()

			// define restore function
			restore := func() {
				if it == 0 {
					for k := 0; k < nip; k++ {
						e.U.States[k].Set(u_states[k])
						e.P.States[k].Set(p_states[k])
					}
					return
				}
				for k := 0; k < nip; k++ {
					e.U.States[k].Set(u_statesBkp[k])
					e.P.States[k].Set(p_statesBkp[k])
				}
			}

			// check
			o.check("Kuu", d, e, e.U.Umap, e.U.Umap, e.U.K, restore)
			o.check("Kup", d, e, e.U.Umap, e.P.Pmap, e.Kup, restore)
			o.check("Kpu", d, e, e.P.Pmap, e.U.Umap, e.Kpu, restore)
			o.check("Kpp", d, e, e.P.Pmap, e.P.Pmap, e.P.Kpp, restore)
			o.check("Kpf", d, e, e.P.Pmap, e.P.Fmap, e.P.Kpf, restore)
			o.check("Kfp", d, e, e.P.Fmap, e.P.Pmap, e.P.Kfp, restore)
			o.check("Kff", d, e, e.P.Fmap, e.P.Fmap, e.P.Kff, restore)
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

			// copy states and solution
			nip := len(e.Rod.IpsElem)
			states := make([]*msolid.OnedState, nip)
			statesBkp := make([]*msolid.OnedState, nip)
			for i := 0; i < nip; i++ {
				states[i] = e.States[i].GetCopy()
				statesBkp[i] = e.StatesBkp[i].GetCopy()
			}
			o.aux_arrays(d)

			// make sure to restore states and solution
			defer func() {
				for i := 0; i < nip; i++ {
					e.States[i].Set(states[i])
					e.StatesBkp[i].Set(statesBkp[i])
				}
				copy(d.Sol.ΔY, o.ΔYbkp)
			}()

			// define restore function
			restore := func() {
				if it == 0 {
					for k := 0; k < nip; k++ {
						e.States[k].Set(states[k])
					}
					return
				}
				for k := 0; k < nip; k++ {
					e.States[k].Set(statesBkp[k])
				}
			}

			// check
			o.check("Krr", d, e, e.Rod.Umap, e.Rod.Umap, e.Krr, restore)
			o.check("Krs", d, e, e.Rod.Umap, e.Sld.Umap, e.Krs, restore)
			o.check("Ksr", d, e, e.Sld.Umap, e.Rod.Umap, e.Ksr, restore)
			o.check("Kss", d, e, e.Sld.Umap, e.Sld.Umap, e.Kss, restore)
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

			// copy states and solution
			nip := len(e.LinIps)
			states := make([]*msolid.OnedState, nip)
			statesBkp := make([]*msolid.OnedState, nip)
			for i := 0; i < nip; i++ {
				states[i] = e.States[i].GetCopy()
				statesBkp[i] = e.StatesBkp[i].GetCopy()
			}
			o.aux_arrays(d)

			// make sure to restore states and solution
			defer func() {
				for i := 0; i < nip; i++ {
					e.States[i].Set(states[i])
					e.StatesBkp[i].Set(statesBkp[i])
				}
				copy(d.Sol.ΔY, o.ΔYbkp)
			}()

			// define restore function
			restore := func() {
				if it == 0 {
					for k := 0; k < nip; k++ {
						e.States[k].Set(states[k])
					}
					return
				}
				for k := 0; k < nip; k++ {
					e.States[k].Set(statesBkp[k])
				}
			}

			// check
			o.check("Kll", d, e, e.LinUmap, e.LinUmap, e.Kll, restore)
			o.check("Kls", d, e, e.LinUmap, e.SldUmap, e.Kls, restore)
			o.check("Ksl", d, e, e.SldUmap, e.LinUmap, e.Ksl, restore)
			o.check("Kss", d, e, e.SldUmap, e.SldUmap, e.Kss, restore)
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
func (o *testKb) aux_arrays(d *Domain) {
	if len(o.Fbtmp) != d.Ny {
		o.Fbtmp = make([]float64, d.Ny)
		o.Yold = make([]float64, d.Ny)
		o.ΔYbkp = make([]float64, d.Ny)
	}
	for i := 0; i < d.Ny; i++ {
		o.Yold[i] = d.Sol.Y[i] - d.Sol.ΔY[i]
		o.ΔYbkp[i] = d.Sol.ΔY[i]
	}
}

// check performs the checking of Kb using numerical derivatives
func (o *testKb) check(label string, d *Domain, e Elem, Imap, Jmap []int, Kana [][]float64, restore func()) {
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
	step := 1e-6
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
				restore()
				if ele, ok := e.(ElemIntvars); ok {
					ele.Update(d.Sol)
				}
				e.AddToRhs(o.Fbtmp, d.Sol)
				res = -o.Fbtmp[I]
				d.Sol.Y[J] = tmp
				return res
			}, d.Sol.Y[J], step)
			chk.AnaNum(o.tst, io.Sf(label+"%3d%3d", i, j), o.tol, Kana[i][j], dnum, o.verb)
		}
	}
}
