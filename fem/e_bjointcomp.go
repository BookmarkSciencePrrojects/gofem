// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/msolid"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
)

// BjointComp implements a beam-joint (interface/link) element for embedded beams with nodes
// compatible with the nodes of the surrounding solid elements
//  Note: beamNu corresponds to the number of displacemetns DOFs of beam; i.e. without rotations
type BjointComp struct {

	// basic data
	Sim    *inp.Simulation // simulation
	Cell   *inp.Cell       // the cell structure
	Edat   *inp.ElemData   // element data; stored in allocator to be used in Connect
	Ny     int             // total number of dofs == LinNu + sld.Nu where LinNu is the number of displacement DOFs of beam == 2*Ndim
	Ndim   int             // space dimension
	TolNod float64         // tolerance to find beam/solid compatible nodes

	// essential
	Lin *Beam           // beam (line) element
	Sld *ElemU          // solid element
	Mdl msolid.RjointM1 // material model

	// additional
	SldLocVid []int // local vertices IDs of solid nodes connected to beam through this joint

	// shape and integration points
	LinShp *shp.Shape   // lin2 shape
	LinIps []shp.Ipoint // integration points along line of beam / joint

	// parameters
	h  float64 // perimeter of beam element
	k1 float64 // lateral stiffness
	k2 float64 // lateral stiffness

	// variables for Coulomb model
	Coulomb bool      // use Coulomb model
	t1      []float64 // [ndim] traction vectors for σc
	t2      []float64 // [ndim] traction vectors for σc

	// auxiliary variables
	fC []float64 // [beamNu] internal/contact forces vector
	qb []float64 // [ndim] resultant traction vector 'holding' the beam @ ip
	/*
		ΔuC [][]float64 // [rodNn][ndim] relative displ. increment of solid @ nodes of beam
		Δw  []float64   // [ndim] relative velocity
	*/

	// temporary Jacobian matrices
	Kll [][]float64 // [linNu][linNu] K_lin_lin
	Kls [][]float64 // [linNu][sldNu] K_lin_sld
	Ksl [][]float64 // [sldNu][linNu] K_sld_lin
	Kss [][]float64 // [sldNu][sldNu] K_sld_sld

	// internal values
	States    []*msolid.OnedState // [nip] internal states
	StatesBkp []*msolid.OnedState // [nip] backup internal states
	StatesAux []*msolid.OnedState // [nip] backup internal states
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	infogetters["bjointcomp"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData) *Info {
		return &Info{}
	}

	// element allocator
	eallocators["bjointcomp"] = func(sim *inp.Simulation, cell *inp.Cell, edat *inp.ElemData, x [][]float64) Elem {
		var o BjointComp
		o.Sim = sim
		o.Cell = cell
		o.Edat = edat
		o.Ndim = sim.Ndim
		o.TolNod = 1e-7
		o.LinShp = shp.Get("lin2", cell.GoroutineId)
		if o.LinShp == nil {
			chk.Panic("cannot allocate \"lin2\" shape for beam/joint {tag=%d id=%d material=%q}", cell.Tag, cell.Id, edat.Mat)
		}
		var err error
		o.LinIps, _, err = o.LinShp.GetIps(edat.Nip, 0)
		if err != nil {
			chk.Panic("cannot get integration points for beam/joint {tag=%d id=%d material=%q} with nip=%d", cell.Tag, cell.Id, edat.Mat, edat.Nip)
		}
		return &o
	}
}

// Id returns the cell Id
func (o *BjointComp) Id() int { return o.Cell.Id }

// Connect connects rod/solid elements in this BjointComp
func (o *BjointComp) Connect(cid2elem []Elem, cell *inp.Cell) (nnzK int, err error) {

	// get beam and solid elements
	linId := cell.JlinId
	sldId := cell.JsldId
	o.Lin = cid2elem[linId].(*Beam)
	o.Sld = cid2elem[sldId].(*ElemU)
	if o.Lin == nil {
		err = chk.Err("cannot find joint's beam cell with id == %d", linId)
		return
	}
	if o.Sld == nil {
		err = chk.Err("cannot find joint's solid cell with id == %d", sldId)
		return
	}

	// number of nodes of beam/line and DOFs
	linNn := 2
	linNu := linNn * o.Ndim // displacement DOFs only
	o.Ny = linNu + o.Sld.Nu // total number of DOFs

	// data from solid element
	sldH := o.Sld.Cell.Shp
	sldNn := sldH.Nverts
	sldNu := o.Sld.Nu

	// local vertex ids of solid vertices connected to beam via joint
	o.SldLocVid = make([]int, 2)
	var idx int
	var dist float64
	for i := 0; i < linNn; i++ {
		for j := 0; j < sldNn; j++ {
			dist = 0.0
			for k := 0; k < o.Ndim; k++ {
				dist += math.Pow(o.Lin.X[k][i]-o.Sld.X[k][j], 2.0)
			}
			dist = math.Sqrt(dist)
			if dist < o.TolNod {
				o.SldLocVid[idx] = j
				idx++
				break
			}
		}
	}
	io.Pforan("SldLocVid = %v\n", o.SldLocVid)

	// material model name
	matdata := o.Sim.MatParams.Get(o.Edat.Mat)
	if matdata == nil {
		err = chk.Err("materials database failed on getting %q material\n", o.Edat.Mat)
		return
	}

	// initialise model
	err = o.Mdl.Init(matdata.Prms)
	if err != nil {
		err = chk.Err("model initialisation failed:\n%v", err)
		return
	}

	// parameters
	for _, p := range matdata.Prms {
		switch p.N {
		case "h":
			o.h = p.V
		case "k1":
			o.k1 = p.V
		case "k2":
			o.k2 = p.V
		case "mu": // no need to store value here because it goes to model
			if p.V > 0.0 {
				o.Coulomb = true
			}
		}
	}

	// auxiliary variables
	o.fC = make([]float64, linNu)
	o.qb = make([]float64, o.Ndim)
	/*
		o.ΔuC = la.MatAlloc(2, o.Ndim)
		o.Δw = make([]float64, o.Ndim)
	*/

	// temporary Jacobian matrices. see Eq. (57)
	o.Kll = la.MatAlloc(linNu, linNu)
	o.Kls = la.MatAlloc(linNu, sldNu)
	o.Ksl = la.MatAlloc(sldNu, linNu)
	o.Kss = la.MatAlloc(sldNu, sldNu)

	// success
	return o.Ny * o.Ny, nil
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// SetEqs set equations
func (o *BjointComp) SetEqs(eqs [][]int, mixedform_eqs []int) (err error) {
	return
}

// SetEleConds set element conditions
func (o *BjointComp) SetEleConds(key string, f fun.Func, extra string) (err error) {
	return
}

// InterpStarVars interpolates star variables to integration points
func (o *BjointComp) InterpStarVars(sol *Solution) (err error) {
	return
}

// adds -R to global residual vector fb
func (o *BjointComp) AddToRhs(fb []float64, sol *Solution) (err error) {

	// internal forces vector
	la.VecFill(o.fC, 0)

	// auxiliary
	linNn := 2                  // number of nodes
	linNdof := 3 * (o.Ndim - 1) // number of DOFs, including rotational ones
	e0, e1, e2 := o.Lin.e0, o.Lin.e1, o.Lin.e2

	// loop over integration points along line
	var coef, τ, qn1, qn2 float64
	for idx, ip := range o.LinIps {

		// interpolation functions and gradients
		err = o.LinShp.CalcAtIp(o.Lin.X, ip, true)
		if err != nil {
			return
		}
		coef = ip[3] * o.LinShp.J
		S := o.LinShp.S

		// state variables
		τ = o.States[idx].Sig
		qn1 = o.States[idx].Phi[0]
		qn2 = o.States[idx].Phi[1]

		// fC vector
		for i := 0; i < o.Ndim; i++ {
			o.qb[i] = τ*o.h*e0[i] + qn1*e1[i] + qn2*e2[i]
			for m := 0; m < linNn; m++ {
				r := i + m*o.Ndim
				o.fC[r] += coef * S[m] * o.qb[i]
			}
		}
	}

	// fb = -Resid;  fB = -fC  and  fS = fC  =>  fb := {fC, -fC}
	var r, s, I int
	for i := 0; i < o.Ndim; i++ {
		for m := 0; m < linNn; m++ {
			r = i + m*o.Ndim  // index in fC: does not consider rotation DOFs
			s = i + m*linNdof // index in beam's Umap: considers rotations
			I = o.Lin.Umap[s] // global index
			fb[I] += o.fC[r]  // fb := -fB = +fC
		}
		for _, m := range o.SldLocVid {
			r = i + m*o.Ndim  // index in fC: no rotations
			I = o.Sld.Umap[r] // global index
			fb[I] -= o.fC[r]  // fb := -fS = -fC
		}
	}
	return
}

// adds element K to global Jacobian matrix Kb
func (o *BjointComp) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (err error) {

	// auxiliary
	linNn := 2
	linNu := linNn * o.Ndim
	sldH := o.Sld.Cell.Shp
	sldNn := sldH.Nverts
	e0, e1, e2 := o.Lin.e0, o.Lin.e1, o.Lin.e2

	// zero K matrices
	for i := 0; i < linNu; i++ {
		for j := 0; j < linNu; j++ {
			o.Kll[i][j] = 0
		}
		for j := 0; j < o.Sld.Nu; j++ {
			o.Kls[i][j] = 0
			o.Ksl[j][i] = 0
		}
	}
	la.MatFill(o.Kss, 0)

	// auxiliary
	var coef float64
	var DτDω float64
	var Dwb0Du_nj, Dwb1Du_nj, Dwb2Du_nj float64
	var DτDu_nj, DqbDu_nij float64
	var Dwb0Dur_nj, Dwb1Dur_nj, Dwb2Dur_nj float64
	var DqbDur_nij float64

	// loop over line's integration points
	for idx, ip := range o.LinIps {

		// interpolation functions and gradients
		err = o.LinShp.CalcAtIp(o.Lin.X, ip, true)
		if err != nil {
			return
		}
		coef = ip[3] * o.LinShp.J
		S := o.LinShp.S

		// model derivatives
		DτDω, err = o.Mdl.CalcD(o.States[idx], firstIt)
		if err != nil {
			return
		}

		// compute derivatives
		for j := 0; j < o.Ndim; j++ {

			// Krr and Ksl; derivatives with respect to ur_nj
			for n := 0; n < linNn; n++ {

				// ∂wb/∂ur Eq (A.4)
				Dwb0Dur_nj = -S[n] * e0[j]
				Dwb1Dur_nj = -S[n] * e1[j]
				Dwb2Dur_nj = -S[n] * e2[j]

				// compute ∂■/∂ur derivatives
				c := j + n*o.Ndim
				for i := 0; i < o.Ndim; i++ {

					// ∂qb/∂ur
					DqbDur_nij = o.h*e0[i]*(DτDω*Dwb0Dur_nj) + o.k1*e1[i]*Dwb1Dur_nj + o.k2*e2[i]*Dwb2Dur_nj

					// Krr := ∂fr/∂ur
					for m := 0; m < linNn; m++ {
						r := i + m*o.Ndim
						o.Kll[r][c] -= coef * S[m] * DqbDur_nij
					}

					//  Ksl := ∂fs/∂ur
					for m := 0; m < linNn; m++ {
						r := i + m*o.Ndim
						o.Ksl[r][c] += coef * S[m] * DqbDur_nij
					}
				}
			}

			// Kls and Kss
			for n := 0; n < sldNn; n++ {

				// ∂wb/∂us
				Dwb0Du_nj, Dwb1Du_nj, Dwb2Du_nj = 0, 0, 0
				for m := 0; m < linNn; m++ {
					Dwb0Du_nj += S[m] * e0[j]
					Dwb1Du_nj += S[m] * e1[j]
					Dwb2Du_nj += S[m] * e2[j]
				}

				// ∂τ/∂us_nj
				DτDu_nj = DτDω * Dwb0Du_nj

				// compute ∂■/∂us derivatives
				c := j + n*o.Ndim
				for i := 0; i < o.Ndim; i++ {

					// ∂qb/∂us
					DqbDu_nij = o.h*e0[i]*DτDu_nj + o.k1*e1[i]*Dwb1Du_nj + o.k2*e2[i]*Dwb2Du_nj

					// Kls := ∂fr/∂us
					for m := 0; m < linNn; m++ {
						r := i + m*o.Ndim
						o.Kls[r][c] -= coef * S[m] * DqbDu_nij
					}

					// Kss := ∂fs/∂us
					for m := 0; m < sldNn; m++ {
						r := i + m*o.Ndim
						o.Kss[r][c] += coef * S[m] * DqbDu_nij
					}
				}
			}
		}
	}

	// add K to sparse matrix Kb
	/*
		for i, I := range o.Rod.Umap {
			for j, J := range o.Rod.Umap {
				Kb.Put(I, J, o.Krr[i][j])
			}
			for j, J := range o.Sld.Umap {
				Kb.Put(I, J, o.Kls[i][j])
				Kb.Put(J, I, o.Ksl[j][i])
			}
		}
		for i, I := range o.Sld.Umap {
			for j, J := range o.Sld.Umap {
				Kb.Put(I, J, o.Kss[i][j])
			}
		}
	*/
	return
}

// Update perform (tangent) update
func (o *BjointComp) Update(sol *Solution) (err error) {
	return
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// Ipoints returns the real coordinates of integration points [nip][ndim]
func (o *BjointComp) Ipoints() (coords [][]float64) {
	coords = la.MatAlloc(len(o.LinIps), o.Ndim)
	for idx, ip := range o.LinIps {
		coords[idx] = o.LinShp.IpRealCoords(o.Lin.X, ip)
	}
	return
}

// SetIniIvs sets initial ivs for given values in sol and ivs map
func (o *BjointComp) SetIniIvs(sol *Solution, ivs map[string][]float64) (err error) {
	nip := len(o.LinIps)
	o.States = make([]*msolid.OnedState, nip)
	o.StatesBkp = make([]*msolid.OnedState, nip)
	o.StatesAux = make([]*msolid.OnedState, nip)
	for i := 0; i < nip; i++ {
		o.States[i], _ = o.Mdl.InitIntVars()
		o.StatesBkp[i] = o.States[i].GetCopy()
		o.StatesAux[i] = o.States[i].GetCopy()
	}
	return
}

// BackupIvs create copy of internal variables
func (o *BjointComp) BackupIvs(aux bool) (err error) {
	if aux {
		for i, s := range o.StatesAux {
			s.Set(o.States[i])
		}
		return
	}
	for i, s := range o.StatesBkp {
		s.Set(o.States[i])
	}
	return
}

// RestoreIvs restore internal variables from copies
func (o *BjointComp) RestoreIvs(aux bool) (err error) {
	if aux {
		for i, s := range o.States {
			s.Set(o.StatesAux[i])
		}
		return
	}
	for i, s := range o.States {
		s.Set(o.StatesBkp[i])
	}
	return
}

// Ureset fixes internal variables after u (displacements) have been zeroed
func (o *BjointComp) Ureset(sol *Solution) (err error) {
	return
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o *BjointComp) Encode(enc Encoder) (err error) {
	return enc.Encode(o.States)
}

// Decode decodes internal variables
func (o *BjointComp) Decode(dec Decoder) (err error) {
	err = dec.Decode(&o.States)
	if err != nil {
		return
	}
	return o.BackupIvs(false)
}

// OutIpsData returns data from all integration points for output
func (o *BjointComp) OutIpsData() (data []*OutIpData) {
	for idx, ip := range o.LinIps {
		s := o.States[idx]
		x := o.LinShp.IpRealCoords(o.Lin.X, ip)
		calc := func(sol *Solution) (vals map[string]float64) {
			vals = make(map[string]float64)
			vals["tau"] = s.Sig
			vals["ompb"] = s.Alp[0]
			return
		}
		data = append(data, &OutIpData{o.Id(), x, calc})
	}
	return
}
