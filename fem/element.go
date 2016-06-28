// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
)

// IpsMap defines a map to hold integration points' results
type IpsMap map[string][]float64

// NewIpsMap returns a new IpsMap
func NewIpsMap() *IpsMap {
	var M IpsMap
	M = make(map[string][]float64)
	return &M
}

// Set sets item in map by key and ip-index. The slice is resized with nip in case it's empty
//  Input:
//   idx -- index of integration point
//   nip -- number of integration points (to resize if necessary)
//   val -- value of 'key' @ integration point 'idx'
func (o *IpsMap) Set(key string, idx, nip int, val float64) {
	if slice, ok := (*o)[key]; ok {
		slice[idx] = val
		return
	}
	slice := make([]float64, nip)
	slice[idx] = val
	(*o)[key] = slice
}

// Get returns item corresponding to 'key' and integration point 'idx'
//  Note: this function returns 0 if 'key' is not found. It also does not check for out-of-bound errors
func (o *IpsMap) Get(key string, idx int) float64 {
	if slice, ok := (*o)[key]; ok {
		return slice[idx]
	}
	return 0
}

// Elem defines what all elements must compute
type Elem interface {

	// information and initialisation
	Id() int                                             // returns the cell Id
	SetEqs(eqs [][]int, mixedform_eqs []int) (err error) // set equations

	// conditions (natural BCs and element's)
	SetEleConds(key string, f fun.Func, extra string) (err error) // set element conditions

	// called for each time step
	InterpStarVars(sol *Solution) (err error) // interpolate star variables to integration points

	// called for each iteration
	AddToRhs(fb []float64, sol *Solution) (err error)                // adds -R to global residual vector fb
	AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (err error) // adds element K to global Jacobian matrix Kb

	// reading and writing of element data
	Encode(enc Encoder) (err error) // encodes internal variables
	Decode(dec Decoder) (err error) // decodes internal variables
}

// ElemIntvars defines elements with {z,q} internal variables
type ElemIntvars interface {
	Update(sol *Solution) (err error)                              // perform (tangent) update
	SetIniIvs(sol *Solution, ivs map[string][]float64) (err error) // sets initial ivs for given values in sol and ivs map
	BackupIvs(aux bool) (err error)                                // create copy of internal variables
	RestoreIvs(aux bool) (err error)                               // restore internal variables from copies
	Ureset(sol *Solution) (err error)                              // fixes internal variables after u (displacements) have been zeroed
}

// ElemConnector defines connector elements; elements that depend upon others
type ElemConnector interface {
	Id() int                                                    // returns the cell Id
	Connect(cid2elem []Elem, c *inp.Cell) (nnzK int, err error) // connect multiple elements; e.g.: connect rod/solid elements in Rjoints
}

// ElemExtrap defines elements with functions to extrapolate internal values
type ElemExtrap interface {
	AddToExt(sol *Solution) (err error) // adds extrapolated values to global array
}

// ElemOutIps defines elements that can output integration points' values
type ElemOutIps interface {
	Id() int                            // returns the cell Id
	OutIpCoords() [][]float64           // coordinates of integration points
	OutIpKeys() []string                // integration points' keys; e.g. "pl", "sl"
	OutIpVals(M *IpsMap, sol *Solution) // integration points' values corresponding to keys
}

// ElemFixedKM defines elements with fixed K,M matrices; to be recomputed if prms are changed
type ElemFixedKM interface {
	Recompute(withM bool) // recompute K and M
}

// Info holds all information required to set a simulation stage
type Info struct {

	// essential
	Dofs [][]string        // solution variables PER NODE. ex for 2 nodes: [["ux", "uy", "rz"], ["ux", "uy", "rz"]]
	Y2F  map[string]string // maps "y" keys to "f" keys. ex: "ux" => "fx", "pl" => "ql"

	// internal Dofs; e.g. for mixed formulations
	NintDofs int // number of internal dofs

	// t1 and t2 variables (time-derivatives of first and second order)
	T1vars []string // "pl"
	T2vars []string // "ux", "uy"

	// required to be extrapolated; e.g. by beam-joints
	Nextrap int // e.g. "nsig"
}

// GetElemInfo returns information about elements/formulations
//  cellType -- e.g. "qua8"
//  elemType -- e.g. "u"
func GetElemInfo(cell *inp.Cell, reg *inp.Region, sim *inp.Simulation) (info *Info, inactive bool, err error) {
	edat := reg.Etag2data(cell.Tag)
	if edat == nil {
		err = chk.Err("cannot get data for element {tag=%d, id=%d}", cell.Tag, cell.Id)
		return
	}
	inactive = edat.Inact
	infogetter, ok := infogetters[edat.Type]
	if !ok {
		err = chk.Err("cannot get info for element {type=%q, tag=%d, id=%d}", edat.Type, cell.Tag, cell.Id)
		return
	}
	info = infogetter(sim, cell, edat)
	if info == nil {
		err = chk.Err("info for element {type=%q, tag=%d, id=%d} is not available", edat.Type, cell.Tag, cell.Id)
	}
	return
}

// NewElem returns a new element from its type; e.g. "p", "u" or "up"
func NewElem(cell *inp.Cell, reg *inp.Region, sim *inp.Simulation) (ele Elem, err error) {
	edat := reg.Etag2data(cell.Tag)
	if edat == nil {
		err = chk.Err("cannot get data for element {tag=%d, id=%d}", cell.Tag, cell.Id)
		return
	}
	allocator, ok := eallocators[edat.Type]
	if !ok {
		err = chk.Err("cannot get allocator for element {type=%q, tag=%d, id=%d}", edat.Type, cell.Tag, cell.Id)
		return
	}
	x := BuildCoordsMatrix(cell, reg.Msh)
	ele = allocator(sim, cell, edat, x)
	if ele == nil {
		err = chk.Err("element {type=%q, tag=%d, id=%d} is not available", edat.Type, cell.Tag, cell.Id)
	}
	return
}

// BuildCoordsMatrix returns the coordinate matrix of a particular Cell
func BuildCoordsMatrix(cell *inp.Cell, msh *inp.Mesh) (x [][]float64) {
	x = la.MatAlloc(msh.Ndim, len(cell.Verts))
	for i := 0; i < msh.Ndim; i++ {
		for j, v := range cell.Verts {
			x[i][j] = msh.Verts[v].C[i]
		}
	}
	return
}

// infogetters holds all available formulations/info; elemType => infogetter
var infogetters = make(map[string]func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *Info)

// eallocators holds all available elements; elemType => eallocator
var eallocators = make(map[string]func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) Elem)
