// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"
	"sort"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mdl/fld"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/utl"
)

// ColLayer holds information of one porous column layer. It computes pressures and densities
// based on the following expressions (maximum liquid saturation and minimum gas saturation)
//
//    ρL  = ρL0 + Cl・(pl - pl0)   thus   dρL/dpl = Cl
//    ρG  = ρG0 + Cg・(pg - pg0)   thus   dρG/dpg = Cg
//    sl  = slmax
//    sg  = 1 - slmax
//    ρ   = nf・sl・ρL  +  nf・sg・ρG  +  (1 - nf)・ρS
//    ns  = 1 - nf
//    σV  = σV0 + ρ・g・(H - z)
//
type ColLayer struct {

	// geometry
	Tags  []int   // tags of cells within this layer
	Zmin  float64 // coordinate (elevation) at bottom of layer
	Zmax  float64 // coordinate (elevation) at top of layer
	Elems []Elem  // elements in layer
	Nodes []*Node // nodes in layer

	// parameters: porous medium
	SlMax float64 // maximum liquid saturation; e.g. 1.0
	Nf0   float64 // initial (constant) porosity
	RhoS0 float64 // initial density of solids

	// parameters: total stress analysis
	TotRho    float64 // density for total stress analyses
	TotStress bool    // total stress analysis

	// additional data
	K0   float64 // coefficient to multiply effective vertical stresses and obtain horizontal effective stresses
	SigV float64 // state @ top of layer

	// auxiliary
	Grav float64    // gravity
	Liq  *fld.Model // liquid model
	Gas  *fld.Model // gas model
}

// Calc_ρ computes ρ (mixture density)
func (o ColLayer) Calc_ρ(ρL, ρG float64) (ρ float64) {
	sl := o.SlMax
	sg := 1.0 - sl
	nf := o.Nf0
	ρS := o.RhoS0
	ρ = nf*sl*ρL + nf*sg*ρG + (1.0-nf)*ρS
	return
}

// Calc computes state @ level z
func (o ColLayer) Calc(z float64) (pl, pg, ρL, ρG, ρ, σV float64) {
	if o.TotStress {
		ρ = o.TotRho
	} else {
		if o.Liq == nil {
			chk.Panic("liquid model must be non-nil in order to compute pressure along column")
		}
		pl, ρL = o.Liq.Calc(z)
		if o.Gas != nil {
			pg, ρG = o.Gas.Calc(z)
		}
		ρ = o.Calc_ρ(ρL, ρG)
	}
	Δz := o.Zmax - z
	σV = o.SigV + ρ*o.Grav*Δz
	return
}

// IniSetPorous sets porous medium initial state (including geostatic or hydrostatic)
func (o *Domain) IniSetPorous(stg *inp.Stage) (err error) {

	// check
	dat := stg.IniPorous
	nlayers := len(dat.Layers)
	if nlayers < 1 {
		return chk.Err("layers must be given by specifying the tags corresponding to each layer")
	}
	if len(dat.Nu) != nlayers {
		if len(dat.K0) != nlayers {
			return chk.Err("either Nu or K0 for each layer must be specified")
		}
	}

	// constants
	ZLARGE := o.Sim.MaxElev + 1.0
	ZSMALL := -1.0 // column base is at z=0

	// initialise layers
	var L ColLayers
	L = make([]*ColLayer, nlayers)
	ndim := o.Sim.Ndim
	nodehandled := make(map[int]bool)
	ctaghandled := make(map[int]bool) // required to make sure all elements were initialised
	for i, tags := range dat.Layers {

		// new layer
		L[i] = new(ColLayer)
		L[i].Tags = tags
		L[i].Zmin = ZLARGE
		L[i].Zmax = ZSMALL

		// for each tag of cells in this layer
		for itag, tag := range tags {

			// check tags
			cells := o.Msh.CellTag2cells[tag]
			if len(cells) < 1 {
				return chk.Err("there are no cells with tag = %d", tag)
			}

			// set nodes and elements and find min and max z-coordinates
			for _, c := range cells {

				// set elements in layer
				L[i].Elems = append(L[i].Elems, o.Cid2elem[c.Id])

				// set nodes in layer
				for _, v := range c.Verts {
					if !nodehandled[v] {
						L[i].Nodes = append(L[i].Nodes, o.Vid2node[v])
					}
					L[i].Zmin = utl.Min(L[i].Zmin, o.Msh.Verts[v].C[ndim-1])
					L[i].Zmax = utl.Max(L[i].Zmax, o.Msh.Verts[v].C[ndim-1])
					nodehandled[v] = true
				}
				ctaghandled[c.Tag] = true
			}

			// material data
			edat := o.Reg.Etag2data(tag)
			mat := o.Sim.MatModels.Get(edat.Mat)
			if mat == nil {
				return chk.Err("cannot get material/model data corresponding to element tag %d / matname = %q", tag, edat.Mat)
			}

			// parameters
			if mat.Por != nil {

				// parameters: porous medium
				if itag == 0 {
					L[i].SlMax = mat.Por.Lrm.SlMax()
					L[i].RhoS0 = mat.Por.RhoS0
					L[i].Nf0 = mat.Por.Nf0
				} else {
					if math.Abs(L[i].SlMax-mat.Por.Lrm.SlMax()) > 1e-15 {
						return chk.Err("all cells/tags in the same layer must have the same material parameters. SlMax: %g != %g", L[i].SlMax, mat.Por.Lrm.SlMax())
					}
					if math.Abs(L[i].RhoS0-mat.Por.RhoS0) > 1e-15 {
						return chk.Err("all cells/tags in the same layer must have the same material parameters. RhoS0: %g != %g", L[i].RhoS0, mat.Por.RhoS0)
					}
					if math.Abs(L[i].Nf0-mat.Por.Nf0) > 1e-15 {
						return chk.Err("all cells/tags in the same layer must have the same material parameters. nf0: %g != %g", L[i].Nf0, mat.Por.Nf0)
					}
				}
			} else {

				// parameters: total stress analysis
				if mat.Sld != nil {
					if itag == 0 {
						L[i].TotRho = mat.Sld.GetRho()
						L[i].TotStress = true
					} else {
						if math.Abs(L[i].TotRho-mat.Sld.GetRho()) > 1e-15 {
							return chk.Err("all cells/tags in the same layer must have the same material parameters. Rho: %g != %g", L[i].TotRho, mat.Sld.GetRho())
						}
					}
				} else {
					return chk.Err("material/model data is not available; porous or solid element must be set")
				}
			}
		}

		// Nu or K0
		if len(dat.Nu) == nlayers {
			L[i].K0 = dat.Nu[i] / (1.0 - dat.Nu[i])
		} else {
			L[i].K0 = dat.K0[i]
		}
		if L[i].K0 < 1e-7 {
			return chk.Err("K0 or Nu is incorect: K0=%g, Nu=%g", L[i].K0, dat.Nu)
		}

		// auxiliary
		L[i].Grav = o.Sim.Grav0
		L[i].Liq = o.Sim.LiqMdl
		L[i].Gas = o.Sim.GasMdl
	}

	// make sure all elements tags were handled
	lins_and_joints := make(map[int]bool)
	for tag, cells := range o.Msh.CellTag2cells {
		solids := true
		for _, cell := range cells {
			if !cell.IsSolid {
				lins_and_joints[cell.Id] = true
				solids = false
			}
		}
		if solids && !ctaghandled[tag] {
			return chk.Err("geost: there are cells not included in any layer: ctag=%d", tag)
		}
	}

	// sort layers from top to bottom
	sort.Sort(L)

	// set previous/top states in layers and compute Sol.Y
	for i, lay := range L {

		// vertical stress @ top of layer
		if i == 0 {
			if math.Abs(o.Sim.Data.Surch) > 0 {
				lay.SigV = -o.Sim.Data.Surch
			}
		} else {
			_, _, _, _, _, lay.SigV = L[i-1].Calc(L[i-1].Zmin)
		}

		// set nodes
		for _, nod := range lay.Nodes {
			z := nod.Vert.C[ndim-1]
			pl, pg, _, _, _, _ := lay.Calc(z)
			dof := nod.GetDof("pl")
			if dof != nil {
				o.Sol.Y[dof.Eq] = pl
			}
			dof = nod.GetDof("pg")
			if dof != nil {
				o.Sol.Y[dof.Eq] = pg
			}
		}

		// set elements
		for _, elem := range lay.Elems {
			if ele, okk := elem.(ElemIntvars); okk {

				// build slices
				coords := ele.(ElemOutIps).OutIpCoords()
				nip := len(coords)
				svT := make([]float64, nip) // total vertical stresses
				for i := 0; i < nip; i++ {
					z := coords[i][ndim-1]
					_, _, _, _, _, σV := lay.Calc(z)
					svT[i] = -σV
				}
				var ivs map[string][]float64
				if lay.TotStress {
					sx := make([]float64, nip)
					sy := make([]float64, nip)
					sz := make([]float64, nip)
					for i := 0; i < nip; i++ {
						sv := svT[i]
						sh := lay.K0 * sv
						sx[i], sy[i], sz[i] = sh, sv, sh
						if ndim == 3 {
							sx[i], sy[i], sz[i] = sh, sh, sv
						}
					}
					ivs = map[string][]float64{"sx": sx, "sy": sy, "sz": sz}
				} else {
					ivs = map[string][]float64{"svT": svT, "K0": []float64{lay.K0}}
				}

				// set element's states
				err = ele.SetIniIvs(o.Sol, ivs)
				if err != nil {
					return chk.Err("element's internal values setting failed:\n%v", err)
				}
			}
		}
	}

	// set initial values of lines and joints
	for cid, _ := range lins_and_joints {
		elem := o.Cid2elem[cid]
		if elem != nil {
			if ele, okk := elem.(ElemIntvars); okk {
				err = ele.SetIniIvs(o.Sol, nil)
				if err != nil {
					return chk.Err("cannot set lines or joints internal values failed:\n%v", err)
				}
			}
		}
	}
	return
}

// sorting /////////////////////////////////////////////////////////////////////////////////////////

// ColLayers is a set of Layer
type ColLayers []*ColLayer

// Len the length of Layers
func (o ColLayers) Len() int {
	return len(o)
}

// Swap swaps two Layers
func (o ColLayers) Swap(i, j int) {
	o[i], o[j] = o[j], o[i]
}

// Less compares Layers: sort from top to bottom
func (o ColLayers) Less(i, j int) bool {
	return o[i].Zmin > o[j].Zmin
}
