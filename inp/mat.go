// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"encoding/json"
	"path/filepath"

	"github.com/cpmech/gofem/mdl/cnd"
	"github.com/cpmech/gofem/mdl/fld"
	"github.com/cpmech/gofem/mdl/gen"
	"github.com/cpmech/gofem/mdl/lrm"
	"github.com/cpmech/gofem/mdl/por"
	"github.com/cpmech/gofem/mdl/sld"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
)

// Material holds material data
type Material struct {

	// input
	Name  string   `json:"name"`  // name of material
	Type  string   `json:"type"`  // type of material; e.g. "sld", "cnd", "lrm", "fld", "por", "gen"
	Model string   `json:"model"` // name of model; e.g. "dp", "vm", "elast", etc.
	Deps  []string `json:"deps"`  // dependencies; other material names. e.g. ["water", "dryair", "solid1", "conduct1", "lreten1"]
	Prms  fun.Prms `json:"prms"`  // prms holds all model parameters for this material

	// derived
	Sld sld.Model  // pointer to solid model
	Cnd cnd.Model  // pointer to conductivity model
	Lrm lrm.Model  // pointer to retention model
	Liq *fld.Model // pointer to liquid model
	Gas *fld.Model // pointer to gas model
	Por *por.Model // pointer to porous model
	Gen gen.Model  // pointer to generic model
}

// Mats holds materials
type MatsData []*Material

// MatDb implements a database of materials
type MatDb struct {

	// input
	Functions FuncsData `json:"functions"` // all functions
	Materials MatsData  `json:"materials"` // all materials

	// derived
	SLD map[string]*Material // subset with materials/models: solids
	CND map[string]*Material // subset with materials/models: coductivities
	LRM map[string]*Material // subset with materials/models: retention models
	LIQ map[string]*Material // subset with materials/models: liquids
	GAS map[string]*Material // subset with materials/models: gases
	POR map[string]*Material // subset with materials/models: porous materials
	GEN map[string]*Material // subset with materials/models: generic materials
}

// Clean cleans resources
func (o *MatDb) Clean() {
	for _, mat := range o.Materials {
		if mat.Sld != nil {
			mat.Sld.Clean()
		}
	}
}

// ReadMat reads all materials data from a .mat JSON file
func ReadMat(dir, fn string, ndim int, pstress bool, H, grav float64) (mdb *MatDb, err error) {

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
	mdb.SLD = make(map[string]*Material)
	mdb.CND = make(map[string]*Material)
	mdb.LRM = make(map[string]*Material)
	mdb.LIQ = make(map[string]*Material)
	mdb.GAS = make(map[string]*Material)
	mdb.POR = make(map[string]*Material)
	mdb.GEN = make(map[string]*Material)
	for _, m := range mdb.Materials {
		switch m.Type {
		case "sld":
			mdb.SLD[m.Name] = m
			continue
		case "cnd":
			mdb.CND[m.Name] = m
			continue
		case "lrm":
			mdb.LRM[m.Name] = m
			continue
		case "fld":
			gas := false
			res := m.Prms.Find("gas")
			if res != nil {
				if res.V > 0 {
					gas = true
				}
			}
			if gas {
				mdb.GAS[m.Name] = m
			} else {
				mdb.LIQ[m.Name] = m
			}
			continue
		case "por":
			mdb.POR[m.Name] = m
			continue
		case "gen":
			mdb.GEN[m.Name] = m
			continue
		default:
			err = chk.Err("material type %q is incorrect; options are \"sld\", \"cnd\", \"lrm\", \"fld\", \"por\" and \"gen\"", m.Type)
			return
		}
	}

	// alloc/init: solids
	for _, m := range mdb.SLD {
		m.Sld, err = sld.New(m.Model)
		if err != nil {
			err = chk.Err("cannot allocate solid model %q / material %q\n%v", m.Model, m.Name, err)
			return
		}
		err = m.Sld.Init(ndim, pstress, m.Prms)
		if err != nil {
			err = chk.Err("cannot initialise solid model %q / material %q\n%v", m.Model, m.Name, err)
			return
		}
	}

	// alloc/init: conductivities
	for _, m := range mdb.CND {
		m.Cnd, err = cnd.New(m.Model)
		if err != nil {
			err = chk.Err("cannot allocate conductivity model %q / material %q\n%v", m.Model, m.Name, err)
			return
		}
		err = m.Cnd.Init(m.Prms)
		if err != nil {
			err = chk.Err("cannot initialise conductivity model %q / material %q\n%v", m.Model, m.Name, err)
			return
		}
	}

	// alloc/init: liquid retention models
	for _, m := range mdb.LRM {
		m.Lrm, err = lrm.New(m.Model)
		if err != nil {
			err = chk.Err("cannot allocate liquid retention model %q / material %q\n%v", m.Model, m.Name, err)
			return
		}
		err = m.Lrm.Init(m.Prms)
		if err != nil {
			err = chk.Err("cannot initialise liquid retention model %q / material %q\n%v", m.Model, m.Name, err)
			return
		}
	}

	// alloc/init: liquids
	for _, m := range mdb.LIQ {
		m.Liq = new(fld.Model)
		m.Liq.Init(m.Prms, H, grav)
	}

	// alloc/init: gases
	for _, m := range mdb.GAS {
		m.Gas = new(fld.Model)
		m.Gas.Init(m.Prms, H, grav)
	}

	// init: porous
	for _, m := range mdb.POR {
		for _, name := range m.Deps {
			if mm, ok := mdb.SLD[name]; ok {
				m.Sld = mm.Sld
			}
			if mm, ok := mdb.CND[name]; ok {
				m.Cnd = mm.Cnd
			}
			if mm, ok := mdb.LRM[name]; ok {
				m.Lrm = mm.Lrm
			}
			if mm, ok := mdb.LIQ[name]; ok {
				m.Liq = mm.Liq
			}
			if mm, ok := mdb.GAS[name]; ok {
				m.Gas = mm.Gas
			}
		}
		m.Por = new(por.Model)
		err = m.Por.Init(m.Prms, m.Cnd, m.Lrm, m.Liq, m.Gas, grav)
		if err != nil {
			err = chk.Err("cannot initialise porous model (material %q):\n%v\n", m.Name, err)
			return
		}
	}

	// alloc/init: generic
	for _, m := range mdb.GEN {
		m.Gen, err = gen.New(m.Model)
		if err != nil {
			err = chk.Err("cannot allocate generic model %q / material %q\n%v", m.Model, m.Name, err)
			return
		}
		err = m.Gen.Init(ndim, m.Prms)
		if err != nil {
			err = chk.Err("cannot initialise generic model %q / material %q\n%v", m.Model, m.Name, err)
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

// String prints one function
func (o *Material) String() string {
	fun.G_extraindent = "        "
	return io.Sf("    {\n      \"name\"  : %q,\n      \"type\"  : %q,\n      \"model\" : %q,\n      \"deps\"  : %q,\n      \"prms\"  : [\n%v\n      ]\n    }", o.Name, o.Type, o.Model, o.Deps, o.Prms)
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
