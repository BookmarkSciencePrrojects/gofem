// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"encoding/json"
	"math"
	"path/filepath"
	"strings"

	"github.com/cpmech/gofem/mconduct"
	"github.com/cpmech/gofem/mporous"
	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gofem/msolid"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
)

// Material holds material data
type Material struct {

	// input
	Name  string   `json:"name"`  // name of material
	Type  string   `json:"type"`  // type of material; e.g. "solid", "conduct", "reten", "porous"
	Model string   `json:"model"` // name of model; e.g. "dp", "vm", "elast", etc.
	Extra string   `json:"extra"` // extra information about this material
	Prms  fun.Prms `json:"prms"`  // prms holds all model parameters for this material

	// derived
	Solid   msolid.Model   // pointer to actual solid model
	Conduct mconduct.Model // pointer to actual conductivity model
	Reten   mreten.Model   // pointer to actual retention model
	Porous  *mporous.Model // pointer to actual porous model
}

// Mats holds materials
type MatsData []*Material

// MatDb implements a database of materials
type MatDb struct {

	// input
	Functions FuncsData `json:"functions"` // all functions
	Materials MatsData  `json:"materials"` // all materials

	// derived
	Solids   map[string]*Material // subset with materials/models: solids
	Conducts map[string]*Material // subset with materials/models: coductivities
	Retens   map[string]*Material // subset with materials/models: retention models
	Porous   map[string]*Material // subset with materials/models: porous materials
	Groups   map[string]*Material // subset with materials/models: groups
}

// Clean cleans resources
func (o *MatDb) Clean() {
	for _, mat := range o.Materials {
		if mat.Solid != nil {
			mat.Solid.Clean()
		}
	}
}

// ReadMat reads all materials data from a .mat JSON file
func ReadMat(dir, fn string, ndim int, pstress bool) (mdb *MatDb, err error) {

	// new database
	mdb = new(MatDb)

	// read file
	b, err := io.ReadFile(filepath.Join(dir, fn))
	if err != nil {
		return nil, err
	}

	// decode
	err = json.Unmarshal(b, mdb)
	if err != nil {
		return
	}

	// subsets
	mdb.Solids = make(map[string]*Material)
	mdb.Conducts = make(map[string]*Material)
	mdb.Retens = make(map[string]*Material)
	mdb.Porous = make(map[string]*Material)
	mdb.Groups = make(map[string]*Material)
	for _, m := range mdb.Materials {
		switch m.Type {
		case "solid":
			mdb.Solids[m.Name] = m
			continue
		case "conduct":
			mdb.Conducts[m.Name] = m
			continue
		case "reten":
			mdb.Retens[m.Name] = m
			continue
		case "porous":
			mdb.Porous[m.Name] = m
			continue
		case "group":
			mdb.Groups[m.Name] = m
			continue
		default:
			err = chk.Err("material type %q is incorrect; options are \"solid\", \"conduct\", \"reten\", and \"porous\"", m.Type)
			return
		}
	}

	// alloc/init: solids
	for _, m := range mdb.Solids {
		m.Solid, err = msolid.New(m.Model)
		if err != nil {
			return
		}
		err = m.Solid.Init(ndim, pstress, m.Prms)
		if err != nil {
			return
		}
	}

	// alloc/init: conducts
	for _, m := range mdb.Conducts {
		m.Conduct, err = mconduct.New(m.Model)
		if err != nil {
			return
		}
		err = m.Conduct.Init(m.Prms)
		if err != nil {
			return
		}
	}

	// alloc/init: retens
	for _, m := range mdb.Retens {
		m.Reten, err = mreten.New(m.Model)
		if err != nil {
			return
		}
		err = m.Reten.Init(m.Prms)
		if err != nil {
			return
		}
	}

	// alloc: porous
	for _, m := range mdb.Porous {
		m.Porous = new(mporous.Model)
	}

	// handle groups
	porous2group := make(map[string]*Material)
	for _, m := range mdb.Groups {
		matnames := strings.Fields(m.Extra)
		for _, name := range matnames {
			if mm, ok := mdb.Solids[name]; ok {
				m.Solid = mm.Solid
			}
			if mm, ok := mdb.Conducts[name]; ok {
				m.Conduct = mm.Conduct
			}
			if mm, ok := mdb.Retens[name]; ok {
				m.Reten = mm.Reten
			}
			if mm, ok := mdb.Porous[name]; ok {
				m.Porous = mm.Porous
				porous2group[name] = m
			}
		}
	}

	// init: porous
	for _, m := range mdb.Porous {
		g := porous2group[m.Name]
		if g == nil {
			err = chk.Err("cannot initialise porous model (%q) because it does not belong to any group", m.Name)
			return
		}
		if g.Conduct == nil {
			err = chk.Err("porous material (%q) in group (%q) must have conductivity model", m.Name, g.Name)
			return
		}
		if g.Reten == nil {
			err = chk.Err("porous material (%q) in group (%q) must have liquid retention model", m.Name, g.Name)
			return
		}
		err = m.Porous.Init(m.Prms, g.Conduct, g.Reten)
		if err != nil {
			return
		}
	}
	return
}

// Get returns a material
//  Note: returns nil if not found
func (o MatDb) Get(name string) *Material {
	for _, mat := range o.Materials {
		if mat.Name == name {
			return mat
		}
	}
	return nil
}

// FluidData finds liquid/gas data in any porous material;
// otherwise returns default values of water and dry air
//  TODO: read pl0 and pg0; default is 0
func (o MatDb) FluidData() (RhoL0, RhoG0, pl0, pg0, Cl, Cg float64) {
	RhoL0 = 1.0    // [Mg/m³]
	RhoG0 = 0.0012 // [Mg/m³]
	Cl = 4.53e-7   // [Mg/(m³・kPa)]
	Cg = 1.17e-5   // [Mg/(m³・kPa)]
	found := false
	for _, m := range o.Porous {
		if found == false {
			RhoL0 = m.Porous.RhoL0
			RhoG0 = m.Porous.RhoG0
			Cl = m.Porous.Cl
			Cg = m.Porous.Cg
			found = true
		} else {
			if math.Abs(RhoL0-m.Porous.RhoL0) > 1e-15 {
				chk.Panic("properties of fluids must be the same in all porous materials")
			}
			if math.Abs(RhoG0-m.Porous.RhoG0) > 1e-15 {
				chk.Panic("properties of fluids must be the same in all porous materials")
			}
			if math.Abs(Cl-m.Porous.Cl) > 1e-15 {
				chk.Panic("properties of fluids must be the same in all porous materials")
			}
			if math.Abs(Cg-m.Porous.Cg) > 1e-15 {
				chk.Panic("properties of fluids must be the same in all porous materials")
			}
		}
	}
	return
}

// String prints one function
func (o *Material) String() string {
	fun.G_extraindent = "  "
	fun.G_openbrackets = false
	return io.Sf("    {\n      \"name\"  : %q,\n      \"type\"  : %q,\n      \"model\" : %q,\n      \"extra\" : %q,\n      \"prms\"  : [\n%v\n    }", o.Name, o.Type, o.Model, o.Extra, o.Prms)
}

// String prints materials
func (o MatsData) String() string {
	l := "  \"materials\" : [\n"
	for i, m := range o {
		if i > 0 {
			l += ",\n"
		}
		l += io.Sf("%v", m)
	}
	l += "\n  ]"
	return l
}

// String outputs all materials
func (o MatDb) String() string {
	return io.Sf("{\n%v,\n%v\n}", o.Functions, o.Materials)
}
