// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"encoding/json"
	"path/filepath"

	"github.com/cpmech/gofem/mdl/conduct"
	"github.com/cpmech/gofem/mdl/diffusion"
	"github.com/cpmech/gofem/mdl/fluid"
	"github.com/cpmech/gofem/mdl/generic"
	"github.com/cpmech/gofem/mdl/porous"
	"github.com/cpmech/gofem/mdl/retention"
	"github.com/cpmech/gofem/mdl/solid"
	"github.com/cpmech/gofem/mdl/thermomech"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
)

// Material holds material data
type Material struct {

	// input
	Name  string     `json:"name"`  // name of material
	Type  string     `json:"type"`  // type of material; e.g. "gen", "sld", "fld", "cnd", "lrm", "dif", "trm", "por"
	Model string     `json:"model"` // name of model; e.g. "dp", "vm", "elast", etc.
	Deps  []string   `json:"deps"`  // dependencies; other material names. e.g. ["water", "dryair", "solid1", "conduct1", "lreten1"]
	Prms  fun.Params `json:"prms"`  // prms holds all model parameters for this material

	// derived
	Gen generic.Model    // pointer to generic model
	Sld solid.Model      // pointer to solid model
	Liq *fluid.Model     // pointer to liquid model
	Gas *fluid.Model     // pointer to gas model
	Cnd conduct.Model    // pointer to conductivity model
	Lrm retention.Model  // pointer to retention model
	Dif diffusion.Model  // pointer to diffusion model
	Trm thermomech.Model // pointer to thermo-mechanical model
	Por *porous.Model    // pointer to porous model
}

// Mats holds materials
type MatsData []*Material

// MatDb implements a database of materials
type MatDb struct {

	// input
	Functions FuncsData `json:"functions"` // all functions
	Materials MatsData  `json:"materials"` // all materials

	// derived
	GEN map[string]*Material // subset with materials/models: generic materials
	SLD map[string]*Material // subset with materials/models: solids
	LIQ map[string]*Material // subset with materials/models: liquids
	GAS map[string]*Material // subset with materials/models: gases
	CND map[string]*Material // subset with materials/models: coductivities
	LRM map[string]*Material // subset with materials/models: retention models
	DIF map[string]*Material // subset with materials/models: diffusion
	TRM map[string]*Material // subset with materials/models: thermomech
	POR map[string]*Material // subset with materials/models: porous materials
}

// Free frees memory
func (o *MatDb) Free() {
	for _, mat := range o.Materials {
		if mat.Sld != nil {
			mat.Sld.Free()
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
	mdb.GEN = make(map[string]*Material)
	mdb.SLD = make(map[string]*Material)
	mdb.LIQ = make(map[string]*Material)
	mdb.GAS = make(map[string]*Material)
	mdb.CND = make(map[string]*Material)
	mdb.LRM = make(map[string]*Material)
	mdb.DIF = make(map[string]*Material)
	mdb.TRM = make(map[string]*Material)
	mdb.POR = make(map[string]*Material)
	for _, m := range mdb.Materials {
		switch m.Type {
		case "gen":
			mdb.GEN[m.Name] = m
			continue
		case "sld":
			mdb.SLD[m.Name] = m
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
		case "cnd":
			mdb.CND[m.Name] = m
			continue
		case "lrm":
			mdb.LRM[m.Name] = m
			continue
		case "dif":
			mdb.DIF[m.Name] = m
			continue
		case "trm":
			mdb.TRM[m.Name] = m
			continue
		case "por":
			mdb.POR[m.Name] = m
			continue
		default:
			err = chk.Err("material type %q is incorrect; options are \"gem\", \"sld\", \"fld\", \"cnd\", \"lrm\", \"dif\", \"trm\" and \"por\"", m.Type)
			return
		}
	}

	// alloc/init: generic
	for _, m := range mdb.GEN {
		m.Gen, err = generic.New(m.Model)
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

	// alloc/init: solids
	for _, m := range mdb.SLD {
		m.Sld, err = solid.New(m.Model)
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

	// alloc/init: liquids
	for _, m := range mdb.LIQ {
		m.Liq = new(fluid.Model)
		m.Liq.Init(m.Prms, H, grav)
	}

	// alloc/init: gases
	for _, m := range mdb.GAS {
		m.Gas = new(fluid.Model)
		m.Gas.Init(m.Prms, H, grav)
	}

	// alloc/init: conductivities
	for _, m := range mdb.CND {
		m.Cnd, err = conduct.New(m.Model)
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
		m.Lrm, err = retention.New(m.Model)
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

	// alloc/init: diffusion models
	for _, m := range mdb.DIF {
		m.Dif, err = diffusion.New(m.Model)
		if err != nil {
			err = chk.Err("cannot allocate diffusion model %q / material %q\n%v", m.Model, m.Name, err)
			return
		}
		err = m.Dif.Init(ndim, m.Prms)
		if err != nil {
			err = chk.Err("cannot initialise diffusion model %q / material %q\n%v", m.Model, m.Name, err)
			return
		}
	}

	// alloc/init: thermo-mechanical models
	for _, m := range mdb.TRM {
		for _, name := range m.Deps {
			if mm, ok := mdb.SLD[name]; ok {
				m.Sld = mm.Sld
			}
		}
		m.Trm, err = thermomech.New(m.Model)
		if err != nil {
			err = chk.Err("cannot allocate thermo-mechanical model %q / material %q\n%v", m.Model, m.Name, err)
			return
		}
		err = m.Trm.Init(ndim, m.Prms)
		if err != nil {
			err = chk.Err("cannot initialise thermo-mechanical model %q / material %q\n%v", m.Model, m.Name, err)
			return
		}
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
		m.Por = new(porous.Model)
		err = m.Por.Init(m.Prms, m.Cnd, m.Lrm, m.Liq, m.Gas, grav)
		if err != nil {
			err = chk.Err("cannot initialise porous model (material %q):\n%v\n", m.Name, err)
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
