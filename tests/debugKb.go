// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package tests

import (
	"testing"

	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/ele/porous"
	"github.com/cpmech/gofem/ele/seepage"
	"github.com/cpmech/gofem/ele/solid"
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/num"
)

// Kb helps on checking Kb matrices
type Kb struct {

	// input (must)
	Tst          *testing.T // testing structure
	Eid          int        // element id
	Tol          float64    // tolerance to compare K's
	Tol2         float64    // another tolerance to compare K's
	Step         float64    // step for finite differences method
	Verb         bool       // verbose: show results
	Ni, Nj       int        // number of i and j components of K to be tested; -1 means all K components
	ItMin, ItMax int        // limits to consider test; -1 means all iterations
	Tmin, Tmax   float64    // limits to consider test; -1 means all times

	// derived
	it    int       // current iteration
	t     float64   // current time
	Fbtmp []float64 // auxiliary array
	Ybkp  []float64 // auxiliary array
	ΔYbkp []float64 // auxiliary array
	Yold  []float64 // auxiliary array
}

// Liquid defines a global function to debug Kb for p-elements
func Liquid(main *fem.Main, o *Kb) {
	main.DebugKb = func(d *fem.Domain, it int) {

		elem := d.Elems[o.Eid]
		if e, ok := elem.(*seepage.Liquid); ok {

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
			o.check("Kpp", d, e, e.Pmap, e.Pmap, e.Kpp, o.Tol)
			o.check("Kpf", d, e, e.Pmap, e.Fmap, e.Kpf, o.Tol)
			o.check("Kfp", d, e, e.Fmap, e.Pmap, e.Kfp, o.Tol)
			o.check("Kff", d, e, e.Fmap, e.Fmap, e.Kff, o.Tol)
		}
	}
}

// LiquidGas defines a global function to debug Kb for pp-elements
func LiquidGas(main *fem.Main, o *Kb) {
	main.DebugKb = func(d *fem.Domain, it int) {

		elem := d.Elems[o.Eid]
		if e, ok := elem.(*seepage.LiquidGas); ok {

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
			o.check("Kll", d, e, e.Plmap, e.Plmap, e.Kll, o.Tol)
			o.check("Klg", d, e, e.Plmap, e.Pgmap, e.Klg, o.Tol)
			o.check("Kgl", d, e, e.Pgmap, e.Plmap, e.Kgl, o.Tol)
			o.check("Kgg", d, e, e.Pgmap, e.Pgmap, e.Kgg, o.Tol)
			o.check("Klf", d, e, e.Plmap, e.Flmap, e.Klf, o.Tol)
			o.check("Kfl", d, e, e.Flmap, e.Plmap, e.Kfl, o.Tol)
			o.check("Kff", d, e, e.Flmap, e.Flmap, e.Kff, o.Tol)
		}
	}
}

// Solid defines a global function to debug Kb for u-elements
func Solid(main *fem.Main, o *Kb) {
	main.DebugKb = func(d *fem.Domain, it int) {

		elem := d.Elems[o.Eid]
		if e, ok := elem.(*solid.Solid); ok {

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
				o.check("Kqq", d, e, e.Qmap, e.Qmap, e.Kqq, o.Tol)
				o.check("Kqu", d, e, e.Qmap, e.Umap, e.Kqu, o.Tol)
				o.check("Kuq", d, e, e.Umap, e.Qmap, e.Kuq, o.Tol)
			}
			o.check("K", d, e, e.Umap, e.Umap, e.K, o.Tol)
		}
	}
	return
}

// SolidLiquid defines a global function to debug Kb for up-elements
func SolidLiquid(main *fem.Main, o *Kb) {
	main.DebugKb = func(d *fem.Domain, it int) {

		elem := d.Elems[o.Eid]
		if e, ok := elem.(*porous.SolidLiquid); ok {

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
			o.check("Kuu", d, e, e.U.Umap, e.U.Umap, e.U.K, o.Tol)
			o.check("Kup", d, e, e.U.Umap, e.P.Pmap, e.Kup, o.Tol2)
			o.check("Kpu", d, e, e.P.Pmap, e.U.Umap, e.Kpu, o.Tol)
			o.check("Kpp", d, e, e.P.Pmap, e.P.Pmap, e.P.Kpp, o.Tol)
			o.check("Kpf", d, e, e.P.Pmap, e.P.Fmap, e.P.Kpf, o.Tol)
			o.check("Kfp", d, e, e.P.Fmap, e.P.Pmap, e.P.Kfp, o.Tol)
			o.check("Kff", d, e, e.P.Fmap, e.P.Fmap, e.P.Kff, o.Tol)
		}
	}
	return
}

// SolidLiquidGas defines a global function to debug Kb for upp-elements
func SolidLiquidGas(main *fem.Main, o *Kb) {
	main.DebugKb = func(d *fem.Domain, it int) {

		elem := d.Elems[o.Eid]
		if e, ok := elem.(*porous.SolidLiquidGas); ok {

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
			o.check("Kll", d, e, e.P.Plmap, e.P.Plmap, e.P.Kll, o.Tol)
			o.check("Klg", d, e, e.P.Plmap, e.P.Pgmap, e.P.Klg, o.Tol)
			o.check("Kgl", d, e, e.P.Pgmap, e.P.Plmap, e.P.Kgl, o.Tol)
			o.check("Kgg", d, e, e.P.Pgmap, e.P.Pgmap, e.P.Kgg, o.Tol)
			o.check("Klf", d, e, e.P.Plmap, e.P.Flmap, e.P.Klf, o.Tol)
			o.check("Kfl", d, e, e.P.Flmap, e.P.Plmap, e.P.Kfl, o.Tol)
			o.check("Kff", d, e, e.P.Flmap, e.P.Flmap, e.P.Kff, o.Tol)

			// {pl,pg,u} versus {pl,pg,u}
			o.check("Kul", d, e, e.U.Umap, e.P.Plmap, e.Kul, o.Tol2)
			o.check("Kug", d, e, e.U.Umap, e.P.Pgmap, e.Kug, o.Tol2)
			o.check("Klu", d, e, e.P.Plmap, e.U.Umap, e.Klu, o.Tol)
			o.check("Kgu", d, e, e.P.Pgmap, e.U.Umap, e.Kgu, o.Tol)

			// {u} versus {u}
			o.check("Kuu", d, e, e.U.Umap, e.U.Umap, e.U.K, o.Tol)
		}
	}
	return
}

// Rjoint defines a global function to debug Kb for rjoint-elements
func Rjoint(main *fem.Main, o *Kb) {
	main.DebugKb = func(d *fem.Domain, it int) {

		elem := d.Elems[o.Eid]
		if e, ok := elem.(*solid.Rjoint); ok {

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
			o.check("Krr", d, e, e.Rod.Umap, e.Rod.Umap, e.Krr, o.Tol)
			o.check("Krs", d, e, e.Rod.Umap, e.Sld.Umap, e.Krs, o.Tol)
			o.check("Ksr", d, e, e.Sld.Umap, e.Rod.Umap, e.Ksr, o.Tol)
			o.check("Kss", d, e, e.Sld.Umap, e.Sld.Umap, e.Kss, o.Tol)
		} else {
			io.Pfred("warning: Eid=%d does not correspond to Rjoint element\n", o.Eid)
		}
	}
	return
}

// Bjointcomp defines a global function to debug Kb for bjoint-compatible elements
func Bjointcomp(main *fem.Main, o *Kb) {
	main.DebugKb = func(d *fem.Domain, it int) {

		elem := d.Elems[o.Eid]
		if e, ok := elem.(*solid.BjointComp); ok {

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
			o.check("Kll", d, e, e.LinUmap, e.LinUmap, e.Kll, o.Tol)
			o.check("Kls", d, e, e.LinUmap, e.SldUmap, e.Kls, o.Tol)
			o.check("Ksl", d, e, e.SldUmap, e.LinUmap, e.Ksl, o.Tol)
			o.check("Kss", d, e, e.SldUmap, e.SldUmap, e.Kss, o.Tol)
		} else {
			io.Pfred("warning: Eid=%d does not correspond to BjointComp element\n", o.Eid)
		}
	}
	return
}

// skip skips test based on it and/or t
func (o *Kb) skip() bool {
	if o.ItMin >= 0 {
		if o.it < o.ItMin {
			return true // skip
		}
	}
	if o.ItMax >= 0 {
		if o.it > o.ItMax {
			return true // skip
		}
	}
	if o.Tmin >= 0 {
		if o.t < o.Tmin {
			return true // skip
		}
	}
	if o.Tmax >= 0 {
		if o.t > o.Tmax {
			return true // skip
		}
	}
	if o.Verb {
		io.PfYel("\nit=%2d t=%v\n", o.it, o.t)
	}
	return false
}

// aux_arrays generates auxiliary arrays
func (o *Kb) aux_backup(d *fem.Domain) {
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

func (o *Kb) aux_restore(d *fem.Domain) {
	for i := 0; i < d.Ny; i++ {
		d.Sol.Y[i] = o.Ybkp[i]
		d.Sol.ΔY[i] = o.ΔYbkp[i]
	}
	for _, eivs := range d.ElemIntvars {
		eivs.RestoreIvs(true) // true => aux
	}
}

// check performs the checking of Kb using numerical derivatives
func (o *Kb) check(label string, d *fem.Domain, e ele.Element, Imap, Jmap []int, Kana [][]float64, tol float64) {
	var imap, jmap []int
	if o.Ni < 0 {
		imap = Imap
	} else {
		if o.Ni <= len(Imap) {
			imap = Imap[:o.Ni]
		}
	}
	if o.Nj < 0 {
		jmap = Jmap
	} else {
		if o.Nj <= len(Jmap) {
			jmap = Jmap[:o.Nj]
		}
	}
	if o.Step < 1e-14 {
		o.Step = 1e-6
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
			}, d.Sol.Y[J], o.Step)
			chk.AnaNum(o.Tst, io.Sf(label+"%3d%3d", i, j), tol, Kana[i][j], dnum, o.Verb)
		}
	}
}
