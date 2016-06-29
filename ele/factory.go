// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package ele

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/chk"
)

// InfoFuncType defines a function that returns information about a certain element type
type InfoFuncType func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *Info

// AllocatorType defines a function that allocates an element
type AllocatorType func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) Element

// GetInfo returns information about elements from factory
func GetInfo(cell *inp.Cell, reg *inp.Region, sim *inp.Simulation) (info *Info, inactive bool, err error) {
	edat := reg.Etag2data(cell.Tag)
	if edat == nil {
		err = chk.Err("cannot get data for element {tag=%d, id=%d}", cell.Tag, cell.Id)
		return
	}
	inactive = edat.Inact
	fcn, ok := infofactory[edat.Type]
	if !ok {
		err = chk.Err("cannot get info for element {type=%q, tag=%d, id=%d}", edat.Type, cell.Tag, cell.Id)
		return
	}
	info = fcn(sim, cell, edat)
	if info == nil {
		err = chk.Err("info for element {type=%q, tag=%d, id=%d} is not available", edat.Type, cell.Tag, cell.Id)
	}
	return
}

// New returns a new element from from factory
func New(cell *inp.Cell, reg *inp.Region, sim *inp.Simulation) (ele Element, err error) {
	edat := reg.Etag2data(cell.Tag)
	if edat == nil {
		err = chk.Err("cannot get data for element {tag=%d, id=%d}", cell.Tag, cell.Id)
		return
	}
	fcn, ok := allocators[edat.Type]
	if !ok {
		err = chk.Err("cannot get allocator for element {type=%q, tag=%d, id=%d}", edat.Type, cell.Tag, cell.Id)
		return
	}
	x := BuildCoordsMatrix(cell, reg.Msh)
	ele = fcn(sim, cell, edat, x)
	if ele == nil {
		err = chk.Err("element {type=%q, tag=%d, id=%d} is not available", edat.Type, cell.Tag, cell.Id)
	}
	return
}

// SetInfoFunc sets a new callback function to return information about an element
func SetInfoFunc(elementName string, fcn InfoFuncType) {
	if _, ok := infofactory[elementName]; ok {
		chk.Panic("cannot set information function for %q because element name exists already", elementName)
	}
	infofactory[elementName] = fcn
}

// SetAllocator sets a new callback function to allocate an element
func SetAllocator(elementName string, fcn AllocatorType) {
	if _, ok := allocators[elementName]; ok {
		chk.Panic("cannot set allocator function for %q because element name exists already", elementName)
	}
	allocators[elementName] = fcn
}

// GetInfoFunc gets callback function to return information about an element
func GetInfoFunc(elementName string) InfoFuncType {
	if fcn, ok := infofactory[elementName]; ok {
		return fcn
	}
	chk.Panic("cannot get function for information about element %q", elementName)
	return nil
}

// GetAllocator gets callback function to allocate an element
func GetAllocator(elementName string) AllocatorType {
	if fcn, ok := allocators[elementName]; ok {
		return fcn
	}
	chk.Panic("cannot get allocator function for element %q", elementName)
	return nil
}

// infofactory holds all functions that return information about an element
var infofactory = make(map[string]InfoFuncType)

// allocators holds all element allocators
var allocators = make(map[string]AllocatorType)
