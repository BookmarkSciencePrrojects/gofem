// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mdl/sld"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/tsr"
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
	Lin *Beam         // beam (line) element
	Sld *ElemU        // solid element
	Mdl *sld.RjointM1 // material model

	// asembly maps
	LinUmap []int // beam umap with displacement DOFs equations only
	SldUmap []int // solid umap with displacement DOFs at nodes connected to beam

	// shape and integration points
	LinShp *shp.Shape   // lin2 shape
	LinIps []shp.Ipoint // integration points along line of beam / joint

	// variables for Coulomb model (scratchpad)
	t1 []float64 // [3] traction vectors for σc
	t2 []float64 // [3] traction vectors for σc
	σ  []float64 // stresses @ ips

	// auxiliary variables
	fC []float64 // [beamNu] internal/contact forces vector
	q  []float64 // [ndim] resultant traction vector 'holding' the beam @ ip
	Δw []float64 // [ndim] relative displacement

	// temporary Jacobian matrices
	Kll [][]float64 // [linNu][linNu] K_lin_lin: ∂fl/∂ub
	Kls [][]float64 // [linNu][linNu] K_lin_sld: ∂fl/∂u
	Ksl [][]float64 // [linNu][linNu] K_sld_lin: ∂fs/∂ub
	Kss [][]float64 // [linNu][linNu] K_sld_sld: ∂fs/∂u

	// internal values
	States    []*sld.OnedState // [nip] internal states
	StatesBkp []*sld.OnedState // [nip] backup internal states
	StatesAux []*sld.OnedState // [nip] backup internal states
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

	// auxiliary
	linNn := 2                  // number of nodes of beam
	linNu := linNn * o.Ndim     // all displacement DOFs only
	nodNdof := 3 * (o.Ndim - 1) // number of DOFs per node, including rotational ones

	// total number of DOFs
	o.Ny = linNu + o.Sld.Nu

	// local vertices IDs of solid nodes connected to beam through this joint
	chk.IntAssert(len(o.Cell.JntConVerts), 2) // global vertices IDs
	sldlocvid := make([]int, linNn)
	for i, a := range o.Cell.JntConVerts {
		for j, b := range o.Sld.Cell.Verts {
			if a == b {
				sldlocvid[i] = j
				continue
			}
		}
	}

	// assembly map with displacements DOFs only of beam nodes
	o.LinUmap = make([]int, linNu)
	for m := 0; m < linNn; m++ {
		for i := 0; i < o.Ndim; i++ {
			r := i + m*o.Ndim
			s := i + m*nodNdof
			o.LinUmap[r] = o.Lin.Umap[s]
		}
	}

	// assembly map with DOFs of solid nodes connected to beam
	o.SldUmap = make([]int, linNu)
	for m, n := range sldlocvid {
		for i := 0; i < o.Ndim; i++ {
			r := i + m*o.Ndim
			s := i + n*o.Ndim
			o.SldUmap[r] = o.Sld.Umap[s]
		}
	}

	// model
	mat := o.Sim.MatModels.Get(o.Edat.Mat)
	if mat == nil {
		err = chk.Err("materials database failed on getting %q material\n", o.Edat.Mat)
		return
	}
	o.Mdl = mat.Sld.(*sld.RjointM1)

	// variables for Coulomb model (scratchpad)
	nsig := 2 * o.Ndim
	o.t1 = make([]float64, 3)
	o.t2 = make([]float64, 3)
	o.σ = make([]float64, nsig)

	// auxiliary variables
	o.fC = make([]float64, linNu)
	o.q = make([]float64, o.Ndim)
	o.Δw = make([]float64, o.Ndim)

	// temporary Jacobian matrices. see Eq. (57)
	o.Kll = la.MatAlloc(linNu, linNu)
	o.Kls = la.MatAlloc(linNu, linNu)
	o.Ksl = la.MatAlloc(linNu, linNu)
	o.Kss = la.MatAlloc(linNu, linNu)

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
	linNn := 2
	h := o.Mdl.A_h
	e0, e1, e2 := o.Lin.e0, o.Lin.e1, o.Lin.e2

	// loop over integration points along line
	var coef, τ, q1, q2 float64
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
		q1 = o.States[idx].Phi[0]
		q2 = o.States[idx].Phi[1]

		// fC vector
		for i := 0; i < o.Ndim; i++ {
			o.q[i] = τ*h*e0[i] + q1*e1[i] + q2*e2[i]
			for m := 0; m < linNn; m++ {
				r := i + m*o.Ndim
				o.fC[r] += coef * S[m] * o.q[i]
			}
		}
	}

	// fb = -Resid;  fB = -fC  and  fS = fC  =>  fb := {fC, -fC}
	for i, I := range o.LinUmap {
		fb[I] += o.fC[i] // fb := -fB = +fC
	}
	for i, I := range o.SldUmap {
		fb[I] -= o.fC[i] // fb := -fS = -fC
	}
	return
}

// adds element K to global Jacobian matrix Kb
func (o *BjointComp) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (err error) {

	// auxiliary
	linNn := 2
	linNu := linNn * o.Ndim
	h := o.Mdl.A_h
	kl := o.Mdl.A_kl
	e0, e1, e2 := o.Lin.e0, o.Lin.e1, o.Lin.e2

	// zero K matrices
	for i := 0; i < linNu; i++ {
		for j := 0; j < linNu; j++ {
			o.Kll[i][j] = 0
			o.Kls[i][j] = 0
			o.Ksl[i][j] = 0
			o.Kss[i][j] = 0
		}
	}

	// auxiliary
	var coef float64
	var DτDω float64
	var Dwb0Du_nj, Dwb1Du_nj, Dwb2Du_nj float64
	var DτDu_nj, DqDu_nij float64
	var Dwb0Dub_nj, Dwb1Dub_nj, Dwb2Dub_nj float64
	var DqDub_nij float64

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
		DτDω, _, err = o.Mdl.CalcD(o.States[idx], firstIt)
		if err != nil {
			return
		}

		// compute derivatives
		for j := 0; j < o.Ndim; j++ {

			// loop of nodes of line/beam
			for n := 0; n < linNn; n++ {

				// ∂wb/∂ub
				Dwb0Dub_nj = -S[n] * e0[j]
				Dwb1Dub_nj = -S[n] * e1[j]
				Dwb2Dub_nj = -S[n] * e2[j]

				// ∂wb/∂u
				Dwb0Du_nj = -Dwb0Dub_nj
				Dwb1Du_nj = -Dwb1Dub_nj
				Dwb2Du_nj = -Dwb2Dub_nj

				// ∂τ/∂u_nj
				DτDu_nj = DτDω * Dwb0Du_nj

				// compute ∂■/∂ub and ∂■/∂u derivatives
				c := j + n*o.Ndim
				for i := 0; i < o.Ndim; i++ {

					// ∂q/∂ub and ∂q/∂u
					DqDub_nij = h*e0[i]*(DτDω*Dwb0Dub_nj) + kl*e1[i]*Dwb1Dub_nj + kl*e2[i]*Dwb2Dub_nj
					DqDu_nij = h*e0[i]*DτDu_nj + kl*e1[i]*Dwb1Du_nj + kl*e2[i]*Dwb2Du_nj

					// K matrices
					for m := 0; m < linNn; m++ {
						r := i + m*o.Ndim
						o.Kll[r][c] -= coef * S[m] * DqDub_nij
						o.Ksl[r][c] += coef * S[m] * DqDub_nij
						o.Kls[r][c] -= coef * S[m] * DqDu_nij
						o.Kss[r][c] += coef * S[m] * DqDu_nij
					}
				}
			}
		}
	}

	// add K to sparse matrix Kb
	for i, I := range o.LinUmap {
		for j, J := range o.LinUmap {
			Kb.Put(I, J, o.Kll[i][j])
		}
		for j, J := range o.SldUmap {
			Kb.Put(I, J, o.Kls[i][j])
			Kb.Put(J, I, o.Ksl[j][i])
		}
	}
	for i, I := range o.SldUmap {
		for j, J := range o.SldUmap {
			Kb.Put(I, J, o.Kss[i][j])
		}
	}
	return
}

// Update perform (tangent) update
func (o *BjointComp) Update(sol *Solution) (err error) {

	// auxiliary
	linNn := 2
	kl := o.Mdl.A_kl
	e0, e1, e2 := o.Lin.e0, o.Lin.e1, o.Lin.e2

	// for each integration point
	var Δwb0, Δwb1, Δwb2, σcb float64
	for idx, ip := range o.LinIps {

		// interpolation functions and gradients
		err = o.LinShp.CalcAtIp(o.Lin.X, ip, true)
		if err != nil {
			return
		}
		S := o.LinShp.S

		// relative displacements @ ip of line/beam
		for i := 0; i < o.Ndim; i++ {
			o.Δw[i] = 0
			for m := 0; m < linNn; m++ {
				r := i + m*o.Ndim
				I := o.LinUmap[r]
				J := o.SldUmap[r]
				o.Δw[i] += S[m] * (sol.ΔY[J] - sol.ΔY[I])
			}
		}

		// relative displacents in the coratational system
		Δwb0, Δwb1, Δwb2 = 0, 0, 0
		for i := 0; i < o.Ndim; i++ {
			Δwb0 += e0[i] * o.Δw[i]
			Δwb1 += e1[i] * o.Δw[i]
			Δwb2 += e2[i] * o.Δw[i]
		}

		// updated confining stress
		σcb, _, _ = o.confining_pressure_ip(sol)

		// update models
		err = o.Mdl.Update(o.States[idx], σcb, Δwb0)
		if err != nil {
			return
		}
		o.States[idx].Phi[0] += kl * Δwb1 // q1
		o.States[idx].Phi[1] += kl * Δwb2 // q2
	}
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
	o.States = make([]*sld.OnedState, nip)
	o.StatesBkp = make([]*sld.OnedState, nip)
	o.StatesAux = make([]*sld.OnedState, nip)
	for i := 0; i < nip; i++ {
		o.States[i], _ = o.Mdl.InitIntVars1D()
		o.StatesBkp[i] = o.States[i].GetCopy()
		o.StatesAux[i] = o.States[i].GetCopy()
	}
	return
}

// confining_pressure_ip computes stresses, tractions and confining pressure @ current ip
// after the shape functions and stresses @ nodes are calculated
func (o *BjointComp) confining_pressure_ip(sol *Solution) (σcb, p1, p2 float64) {

	// current shape function @ ip
	S := o.LinShp.S

	// interpolate stresses from nodes
	nsig := 2 * o.Ndim
	for i := 0; i < nsig; i++ {
		o.σ[i] = 0.0
		for m, vid := range o.Cell.JntConVerts {
			o.σ[i] += S[m] * sol.Ext[vid][i]
		}
	}

	// calculate t1 and t2
	for i := 0; i < 3; i++ {
		o.t1[i], o.t2[i] = 0, 0
		for j := 0; j < 3; j++ {
			o.t1[i] += tsr.M2T(o.σ, i, j) * o.Lin.e1[j]
			o.t2[i] += tsr.M2T(o.σ, i, j) * o.Lin.e2[j]
		}
	}

	// calculate p1 and p2
	for i := 0; i < 3; i++ {
		p1 += o.t1[i] * o.Lin.e1[i]
		p2 += o.t2[i] * o.Lin.e2[i]
	}

	// confining pressure: compressive is positive
	σcb = -(p1 + p2) / 2.0 // no need ramp function because model will do this
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

/*
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
*/
