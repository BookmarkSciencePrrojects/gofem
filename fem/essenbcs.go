// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"
	"sort"
	"strings"

	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/mdl/fluid"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
)

// EssentialBc holds information about essential bounday conditions such as constrained nodes.
// Lagrange multipliers are used to implement both single- and multi-point constraints.
//  In general, essential bcs / constraints are defined by means of:
//
//      A・y = c
//
//  The resulting Kb matrix will then have the following form:
//      _       _
//     |  K  At  | / δy \   / -R - At*λ \
//     |         | |    | = |           |
//     |_ A   0 _| \ δλ /   \  c - A*y  /
//         Kb       δyb          fb
//
type EssentialBc struct {
	Key   string    // key such as 'ux', 'uy', 'rigid', 'incsup', 'hst'
	Eqs   []int     // equations numbers; can be more than one e.g. for inclined support
	ValsA []float64 // values for matrix A
	Fcn   fun.Func  // function that implements the "c" vector in  A・y = c
}

// EbcArray is an array of EssentialBc's
type EbcArray []*EssentialBc

// EssentialBcs implements a structure to record the definition of essential bcs / constraints.
// Each constraint will have a unique Lagrange multiplier index.
type EssentialBcs struct {
	LiqMdl *fluid.Model // for computing hydrostatic conditions
	EqsIni map[int]bool // equations that depend on initial values
	Bcs    EbcArray     // active essential bcs / constraints
	A      la.Triplet   // matrix of coefficients 'A'
	Am     *la.CCMatrix // compressed form of A matrix
}

// Init initialises this structure
func (o *EssentialBcs) Init(liqmdl *fluid.Model) {
	o.LiqMdl = liqmdl
	o.EqsIni = make(map[int]bool)
	o.Bcs = make([]*EssentialBc, 0)
}

// Build builds the structures required for assembling A matrix
//  nλ   -- is the number of essential bcs / constraints == number of Lagrange multipliers
//  nnzA -- is the number of non-zeros in matrix 'A'
func (o *EssentialBcs) Build(ny int) (nλ, nnzA int) {

	// skip if there are no constraints
	nλ = len(o.Bcs)
	if nλ == 0 {
		return
	}

	// sort bcs to make sure all processors will number Lagrange multipliers in the same order
	sort.Sort(o.Bcs)

	// count number of non-zeros in matrix A
	for _, bc := range o.Bcs {
		nnzA += len(bc.ValsA)
	}

	// set matrix A
	o.A.Init(nλ, ny, nnzA)
	for i, bc := range o.Bcs {
		for j, eq := range bc.Eqs {
			o.A.Put(i, eq, bc.ValsA[j])
		}
	}
	o.Am = o.A.ToMatrix(nil)
	return
}

// AddtoRhs adds the essential bcs / constraints terms to the augmented fb vector
func (o *EssentialBcs) AddToRhs(fb []float64, sol *ele.Solution) {

	// skip if there are no constraints
	if len(o.Bcs) == 0 {
		return
	}

	// add -At*λ to fb
	la.SpMatTrVecMulAdd(fb, -1, o.Am, sol.L) // fb += -1 * At * λ

	// assemble -rc = c - A*y into fb
	ny := len(sol.Y)
	for i, bc := range o.Bcs {
		fb[ny+i] = bc.Fcn.F(sol.T, nil)
	}
	la.SpMatVecMulAdd(fb[ny:], -1, o.Am, sol.Y) // fb += -1 * A * y
}

// GetIsEssenKeyMap returns the "YandC" map with special keys that EssentialBcs can handle,
// including:
//  rigid  -- define rigid element constraints
//  incsup -- inclined support constraints
//  hst    -- set hydrostatic pressures
func GetIsEssenKeyMap() map[string]bool {
	return map[string]bool{"rigid": true, "incsup": true, "hst": true}
}

// Set sets a constraint if it does not exist yet.
//  key   -- can be Dof key such as "ux", "uy" or constraint type such as "incsup" or "rigid"
//  extra -- is a keycode-style data. e.g. "!type:incsup2d !alp:30"
//  Notes:
//   1) the default key is single point constraint; e.g. "ux", "uy", ...
//   2) hydraulic head can be set with key == "hst" (hydrostatic). In this case, fcn==shift
//      computes a 'shift' such that:
//          pl(t,z) = pl(z) - shift(t)
//   3) if the key as a suffix "_ini", the initial value of essential key will be multiplied
//      by fcn==mult in order to define the boundary condition according to:
//          y(t,z) = y(z)_ini・mult(t)
func (o *EssentialBcs) Set(key string, nodes []*Node, fcn fun.Func, extra string) (err error) {

	// auxiliary
	chk.IntAssertLessThan(0, len(nodes)) // 0 < len(nod)
	if nodes[0] == nil {
		return
	}
	ndim := len(nodes[0].Vert.C)

	// rigid element
	if key == "rigid" {
		a := nodes[0].Dofs
		for i := 1; i < len(nodes); i++ {
			for j, b := range nodes[i].Dofs {
				o.set_eqs(key, []int{a[j].Eq, b.Eq}, []float64{1, -1}, &fun.Zero)
			}
		}
		return // success
	}

	// inclined support
	if key == "incsup" {

		// check
		if ndim != 2 {
			return chk.Err("inclined support works only in 2D for now")
		}

		// get data
		var α float64
		if val, found := io.Keycode(extra, "alp"); found {
			α = io.Atof(val) * math.Pi / 180.0
		}
		co, si := math.Cos(α), math.Sin(α)

		// set for all nodes
		for _, nod := range nodes {
			eqx := nod.Dofs[0].Eq
			eqy := nod.Dofs[1].Eq
			o.set_eqs(key, []int{eqx, eqy}, []float64{co, si}, &fun.Zero)
		}
		return // success
	}

	// hydraulic head
	if key == "hst" {

		// check
		if o.LiqMdl == nil {
			return chk.Err("cannot apply hydrostatic (hst) boundary condition because liquid model is not available\n")
		}

		// set for all nodes
		for _, nod := range nodes {

			// get DOF
			d := nod.GetDof("pl")
			if d == nil {
				continue // node doesn't have key. ex: pl in qua8/qua4 elements
			}

			// create function
			z := nod.Vert.C[1] // 2D
			if ndim == 3 {
				z = nod.Vert.C[2] // 3D
			}
			plVal, _ := o.LiqMdl.Calc(z)
			pl := fun.Add{
				B: 1, Fb: &fun.Cte{C: plVal},
				A: -1, Fa: fcn,
			}

			// set constraint
			o.set_eqs("pl", []int{d.Eq}, []float64{1}, &pl)
		}
		return // success
	}

	// set with initial value
	if strings.HasSuffix(key, "_ini") {

		// set for all nodes
		kkey := key[:len(key)-4]
		for _, nod := range nodes {

			// get DOF
			d := nod.GetDof(kkey)
			if d == nil {
				continue // node doesn't have key. ex: pl in qua8/qua4 elements
			}

			// create function
			f := fun.Mul{Fa: fcn, Fb: &fun.Cte{}}

			// set constraint
			o.set_eqs(kkey, []int{d.Eq}, []float64{1}, &f)
			o.EqsIni[d.Eq] = true
		}
		return // success
	}

	// single-point constraint
	for _, nod := range nodes {

		// get DOF
		d := nod.GetDof(key)
		if d == nil {
			continue // node doesn't have key. ex: pl in qua8/qua4 elements
		}

		// set constraint
		o.set_eqs(key, []int{d.Eq}, []float64{1}, fcn)
	}

	// success
	return
}

// FixIniVals fixes functions of BCs that depend on initial values
func (o *EssentialBcs) FixIniVals(sol *ele.Solution) {
	for eq, _ := range o.EqsIni {
		for _, bc := range o.Bcs {
			for _, eqOld := range bc.Eqs {
				if eqOld == eq {
					fcn := bc.Fcn.(*fun.Mul)
					f := fcn.Fb.(*fun.Cte)
					f.C = sol.Y[eq]
				}
			}
		}
	}
}

// List returns a simple list logging bcs at time t
func (o *EssentialBcs) List(t float64) (l string) {
	l = "\n==================================================================\n"
	l += io.Sf("%8s%8s%25s%25s\n", "eq", "key", "value @ t=0", io.Sf("value @ t=%g", t))
	l += "------------------------------------------------------------------\n"
	sort.Sort(o.Bcs)
	for _, bc := range o.Bcs {
		l += io.Sf("%8d%8s%25.13f%25.13f\n", bc.Eqs[0], bc.Key, bc.Fcn.F(0, nil), bc.Fcn.F(t, nil))
	}
	l += "==================================================================\n"
	return
}

// auxiliary /////////////////////////////////////////////////////////////////////////////////////////

// set_eqs sets/replace constraint and equations
func (o *EssentialBcs) set_eqs(key string, eqs []int, valsA []float64, fcn fun.Func) {

	// replace existent
	for _, eq := range eqs {
		for _, bc := range o.Bcs {
			for _, eqOld := range bc.Eqs {
				if eqOld == eq {
					bc.Key, bc.Eqs, bc.ValsA, bc.Fcn = key, eqs, valsA, fcn
					return
				}
			}
		}
	}

	// add new
	o.Bcs = append(o.Bcs, &EssentialBc{key, eqs, valsA, fcn})
}

// functions to implement Sort interface
func (o EbcArray) Len() int      { return len(o) }
func (o EbcArray) Swap(i, j int) { o[i], o[j] = o[j], o[i] }
func (o EbcArray) Less(i, j int) bool {
	sort.Ints(o[i].Eqs)
	sort.Ints(o[j].Eqs)
	return o[i].Eqs[0] < o[j].Eqs[0]
}
