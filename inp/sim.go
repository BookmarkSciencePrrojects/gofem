// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package inp implements the input data read from a (.sim) JSON file
package inp

import (
	"encoding/json"
	goio "io"
	"math"
	"os"
	"path/filepath"

	"github.com/cpmech/gofem/mdl/fluid"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/rnd"
	"github.com/cpmech/gosl/utl"
)

// Data holds global data for simulations
type Data struct {

	// global information
	Desc    string `json:"desc"`    // description of simulation
	Matfile string `json:"matfile"` // materials file path
	DirOut  string `json:"dirout"`  // directory for output; e.g. /tmp/gofem
	Encoder string `json:"encoder"` // encoder name; e.g. "gob" "json" "xml"

	// problem definition and options
	Steady    bool    `json:"steady"`    // steady simulation
	Axisym    bool    `json:"axisym"`    // axisymmetric
	Pstress   bool    `json:"pstress"`   // plane-stress
	NoLBB     bool    `json:"nolbb"`     // do not satisfy Ladyženskaja-Babuška-Brezzi condition; i.e. do not use [qua8,qua4] for u-p formulation
	Stat      bool    `json:"stat"`      // activate statistics
	Wlevel    float64 `json:"wlevel"`    // water level; 0 means use max elevation
	Surch     float64 `json:"surch"`     // surcharge load at surface == qn0
	LiqMat    string  `json:"liq"`       // name of liquid material
	GasMat    string  `json:"gas"`       // name of gas material
	ListBcs   bool    `json:"listbcs"`   // list boundary conditions
	WriteSmat bool    `json:"writesmat"` // writes /tmp/gofem_Kb.smat file for debugging global Jacobian matrix. The simulation will be stopped.
}

// LinSolData holds data for linear solvers
type LinSolData struct {
	Name      string `json:"name"`      // "mumps" or "umfpack"
	Symmetric bool   `json:"symmetric"` // use symmetric solver
	Verbose   bool   `json:"verbose"`   // verbose?
	Timing    bool   `json:"timing"`    // show timing statistics
	Ordering  string `json:"ordering"`  // ordering scheme
	Scaling   string `json:"scaling"`   // scaling scheme
}

// SolverData holds FEM solver data
type SolverData struct {

	// nonlinear solver
	Type    string  `json:"type"`    // nonlinear solver type: {imp, exp, rex} => implicit, explicit, Richardson extrapolation
	NmaxIt  int     `json:"nmaxit"`  // number of max iterations
	Atol    float64 `json:"atol"`    // absolute tolerance
	Rtol    float64 `json:"rtol"`    // relative tolerance
	FbTol   float64 `json:"fbtol"`   // tolerance for convergence on fb
	FbMin   float64 `json:"fbmin"`   // minimum value of fb
	DvgCtrl bool    `json:"dvgctrl"` // use divergence control
	NdvgMax int     `json:"ndvgmax"` // max number of continued divergence
	CteTg   bool    `json:"ctetg"`   // use constant tangent (modified Newton) during iterations
	ShowR   bool    `json:"showr"`   // show residual

	// Richardson's extrapolation
	REnogus  bool    `json:"renogus"`  // Richardson extrapolation: no Gustafsson's step control
	REnssmax int     `json:"renssmax"` // Richardson extrapolation: max number of substeps
	REatol   float64 `json:"reatol"`   // Richardson extrapolation: absolute tolerance
	RErtol   float64 `json:"rertol"`   // Richardson extrapolation: relative tolerance
	REmfac   float64 `json:"remfac"`   // Richardson extrapolation: multiplier factor
	REmmin   float64 `json:"remmin"`   // Richardson extrapolation: min multiplier
	REmmax   float64 `json:"remmax"`   // Richardson extrapolation: max multiplier

	// transient analyses
	DtMin      float64 `json:"dtmin"`      // minium value of Dt for transient (θ and Newmark / Dyn coefficients)
	Theta      float64 `json:"theta"`      // θ-method
	ThGalerkin bool    `json:"thgalerkin"` // use θ = 2/3
	ThLiniger  bool    `json:"thliniger"`  // use θ = 0.878

	// dynamics
	Theta1 float64 `json:"theta1"` // Newmark's method parameter
	Theta2 float64 `json:"theta2"` // Newmark's method parameter
	HHT    bool    `json:"hht"`    // use Hilber-Hughes-Taylor method
	HHTalp float64 `json:"hhtalp"` // HHT α parameter

	// combination of coefficients
	ThCombo1 bool `json:"thcombo1"` // use θ=2/3, θ1=5/6 and θ2=8/9 to avoid oscillations

	// constants
	Eps float64 `json:"eps"` // smallest number satisfying 1.0 + ϵ > 1.0

	// derived
	Itol float64 // iterations tolerance
}

// ElemData holds element data
type ElemData struct {

	// input data
	Tag   int    `json:"tag"`   // tag of element
	Mat   string `json:"mat"`   // material name
	Type  string `json:"type"`  // type of element. ex: u, p, up, rod, beam, rjoint
	Nip   int    `json:"nip"`   // number of integration points; 0 => use default
	Nipf  int    `json:"nipf"`  // number of integration points on face; 0 => use default
	Extra string `json:"extra"` // extra flags (in keycode format). ex: "!thick:0.2 !nip:4"
	Inact bool   `json:"inact"` // whether element starts inactive or not

	// auxiliary/internal
	Lbb bool // LBB element
}

// Region holds region data
type Region struct {

	// input data
	Desc      string      `json:"desc"`      // description of region. ex: ground, indenter, etc.
	Mshfile   string      `json:"mshfile"`   // file path of file with mesh data
	ElemsData []*ElemData `json:"elemsdata"` // list of elements data
	AbsPath   bool        `json:"abspath"`   // mesh filename is given in absolute path

	// derived
	Msh *Mesh // the mesh
}

// FaceBc holds face boundary condition
type FaceBc struct {
	Tag   int      `json:"tag"`   // tag of face
	Keys  []string `json:"keys"`  // key indicating type of bcs. ex: qn, pw, ux, uy, uz, wwx, wwy, wwz
	Funcs []string `json:"funcs"` // name of function. ex: zero, load, myfunction1, etc.
	Extra string   `json:"extra"` // extra information. ex: '!λl:10'
}

// SeamBc holds seam (3D edge) boundary condition
type SeamBc struct {
	Tag   int      `json:"tag"`   // tag of seam
	Keys  []string `json:"keys"`  // key indicating type of bcs. ex: qn, pw, ux, uy, uz, wwx, wwy, wwz
	Funcs []string `json:"funcs"` // name of function. ex: zero, load, myfunction1, etc.
	Extra string   `json:"extra"` // extra information. ex: '!λl:10'
}

// NodeBc holds node boundary condition
type NodeBc struct {
	Tag   int      `json:"tag"`   // tag of node
	Keys  []string `json:"keys"`  // key indicating type of bcs. ex: pw, ux, uy, uz, wwx, wwy, wwz
	Funcs []string `json:"funcs"` // name of function. ex: zero, load, myfunction1, etc.
	Extra string   `json:"extra"` // extra information. ex: '!λl:10'
}

// EleCond holds element condition
type EleCond struct {
	Tag   int      `json:"tag"`   // tag of cell/element
	Keys  []string `json:"keys"`  // key indicating type of condition. ex: "g" (gravity), "qn" for beams, etc.
	Funcs []string `json:"funcs"` // name of function. ex: grav, none
	Extra string   `json:"extra"` // extra information. ex: '!λl:10'
}

// TimeControl holds data for defining the simulation time stepping
type TimeControl struct {
	Tf     float64 `json:"tf"`     // final time
	Dt     float64 `json:"dt"`     // time step size (if constant)
	DtOut  float64 `json:"dtout"`  // time step size for output
	DtFcn  string  `json:"dtfcn"`  // time step size (function name)
	DtoFcn string  `json:"dtofcn"` // time step size for output (function name)

	// derived
	DtFunc  fun.Func // time step function
	DtoFunc fun.Func // output time step function
}

// IniPorousData  holds data for setting initial porous media state (e.g. geostatic, hydrostatic)
type IniPorousData struct {
	Nu     []float64 `json:"nu"`     // [nlayers] Poisson's coefficient to compute effective horizontal state for each layer
	K0     []float64 `json:"K0"`     // [nlayers] or Earth pressure coefficient at rest to compute effective horizontal stresses
	Layers [][]int   `json:"layers"` // [nlayers][ntagsInLayer]; e.g. [[-1,-2], [-3,-4]] => 2 layers
}

// IniStressData holds data for setting initial stresses
type IniStressData struct {
	Hom bool    `json:"hom"` // homogeneous stress distribution
	Iso bool    `json:"iso"` // isotropic state
	Psa bool    `json:"psa"` // plane-strain state
	S0  float64 `json:"s0"`  // Iso => stress value to use in homogeneous and isotropic distribution
	Sh  float64 `json:"sh"`  // Psa => horizontal stress
	Sv  float64 `json""sv"`  // Psa => vertical stress
	Nu  float64 `json:"nu"`  // Psa => Poisson's coefficient for plane-strain state
}

// IniFcnData holds data for setting initial solution values such as Y, dYdt and d2Ydt2
type IniFcnData struct {
	File string   `json:"file"` // file with values at each node is given; filename with path is provided
	Fcns []string `json:"fcns"` // functions F(t, x) are given; from functions database
	Dofs []string `json:"dofs"` // degrees of freedom corresponding to "fcns"
}

// IniImportRes holds definitions for importing results from a previous simulation
type IniImportRes struct {
	Dir    string `json:"dir"`    // output directory with previous simulation files
	Fnk    string `json:"fnk"`    // previous simulation file name key (without .sim)
	ResetU bool   `json:"resetu"` // reset/zero u (displacements)
}

// Stage holds stage data
type Stage struct {

	// main
	Desc       string `json:"desc"`       // description of simulation stage. ex: activation of top layer
	Activate   []int  `json:"activate"`   // array of tags of elements to be activated
	Deactivate []int  `json:"deactivate"` // array of tags of elements to be deactivated
	Save       bool   `json:"save"`       // save stage data to binary file
	Load       string `json:"load"`       // load stage data (filename) from binary file
	Skip       bool   `json:"skip"`       // do not run stage

	// specific problems data
	SeepFaces []int          `json:"seepfaces"` // face tags corresponding to seepage faces
	IniPorous *IniPorousData `json:"iniporous"` // initial porous media state (geostatic and hydrostatic included)
	IniStress *IniStressData `json:"inistress"` // initial stress data
	IniFcn    *IniFcnData    `json:"inifcn"`    // set initial solution values such as Y, dYdt and d2Ydt2
	IniImport *IniImportRes  `json:"import"`    // import results from another previous simulation

	// conditions
	EleConds []*EleCond `json:"eleconds"` // element conditions. ex: gravity or beam distributed loads
	FaceBcs  []*FaceBc  `json:"facebcs"`  // face boundary conditions
	SeamBcs  []*SeamBc  `json:"seambcs"`  // seam (3D) boundary conditions
	NodeBcs  []*NodeBc  `json:"nodebcs"`  // node boundary conditions

	// timecontrol
	Control TimeControl `json:"control"` // time control
}

// Simulation holds all simulation data
type Simulation struct {

	// input
	Data      Data       `json:"data"`      // stores global simulation data
	Functions FuncsData  `json:"functions"` // stores all boundary condition functions
	PlotF     *PlotFdata `json:"plotf"`     // plot functions
	Regions   []*Region  `json:"regions"`   // stores all regions
	LinSol    LinSolData `json:"linsol"`    // linear solver data
	Solver    SolverData `json:"solver"`    // FEM solver data
	Stages    []*Stage   `json:"stages"`    // stores all stages

	// derived
	GoroutineId int          // id of goroutine to avoid race problems
	DirOut      string       // directory to save results
	Key         string       // simulation key; e.g. mysim01.sim => mysim01 or mysim01-alias
	EncType     string       // encoder type
	Ndim        int          // space dimension
	MaxElev     float64      // maximum elevation
	Grav0       float64      // gravity constant from stage #0
	MatModels   *MatDb       // materials and models
	LiqMdl      *fluid.Model // liquid model to use when computing density and pressure along column; from stage #0
	GasMdl      *fluid.Model // gas model to use when computing density and pressure along column; from stage #0

	// adjustable parameters
	Adjustable   fun.Prms         // adjustable parameters (not dependent)
	AdjRandom    rnd.Variables    // adjustable parameters that are random variables (not dependent)
	AdjDependent fun.Prms         // adjustable parameters that depend on other adjustable parameters
	adjmap       map[int]*fun.Prm // auxiliary map with adjustable (not dependent)
}

// Simulation //////////////////////////////////////////////////////////////////////////////////////

// Clean cleans resources
func (o *Simulation) Clean() {
	if o.MatModels != nil {
		o.MatModels.Clean()
	}
}

// ReadSim reads all simulation data from a .sim JSON file
func ReadSim(simfilepath, alias string, erasePrev, createDirOut bool, goroutineId int) *Simulation {

	// new sim
	var o Simulation
	o.GoroutineId = goroutineId

	// read file
	b, err := io.ReadFile(simfilepath)
	if err != nil {
		chk.Panic("ReadSim: cannot read simulation file %q", simfilepath)
	}

	// set default values
	o.Solver.SetDefault()
	o.LinSol.SetDefault()

	// decode
	err = json.Unmarshal(b, &o)
	if err != nil {
		chk.Panic("ReadSim: cannot unmarshal simulation file %q", simfilepath)
	}

	// input directory and filename key
	dir := filepath.Dir(simfilepath)
	fn := filepath.Base(simfilepath)
	dir = os.ExpandEnv(dir)
	fnkey := io.FnKey(fn)
	o.Key = fnkey
	if alias != "" {
		o.Key += "-" + alias
	}

	// output directory
	o.DirOut = o.Data.DirOut
	if o.DirOut == "" {
		o.DirOut = "/tmp/gofem/" + fnkey
	}

	// encoder type
	o.EncType = o.Data.Encoder
	if o.EncType != "gob" && o.EncType != "json" {
		o.EncType = "gob"
	}

	// create directory
	if createDirOut {
		err = os.MkdirAll(o.DirOut, 0777)
		if err != nil {
			chk.Panic("cannot create directory for output results (%s): %v", o.DirOut, err)
		}
	}

	// erase previous simulation results
	if erasePrev {
		io.RemoveAll(io.Sf("%s/%s*", o.DirOut, fnkey))
	}

	// set solver constants
	o.Solver.PostProcess()

	// for all regions
	for i, reg := range o.Regions {

		// read mesh
		ddir := dir
		if reg.AbsPath {
			ddir = ""
		}
		reg.Msh, err = ReadMsh(ddir, reg.Mshfile, goroutineId)
		if err != nil {
			chk.Panic("ReadSim: cannot read mesh file:\n%v", err)
		}

		// get ndim and max elevation
		if i == 0 {
			o.Ndim = reg.Msh.Ndim
			o.MaxElev = reg.Msh.MaxElev
		} else {
			if reg.Msh.Ndim != o.Ndim {
				chk.Panic("ReadSim: Ndim value is inconsistent: %d != %d", reg.Msh.Ndim, o.Ndim)
			}
			o.MaxElev = utl.Max(o.MaxElev, reg.Msh.MaxElev)
		}
	}

	// for all stages
	var t float64
	for i, stg := range o.Stages {

		// fix Tf
		if stg.Control.Tf < 1e-14 {
			stg.Control.Tf = 1
		}

		// fix Dt
		if stg.Control.DtFcn == "" {
			if stg.Control.Dt < 1e-14 {
				stg.Control.Dt = 1
			}
			stg.Control.DtFunc = &fun.Cte{C: stg.Control.Dt}
		} else {
			stg.Control.DtFunc, err = o.Functions.Get(stg.Control.DtFcn)
			if err != nil {
				chk.Panic("%v", err)
			}
			stg.Control.Dt = stg.Control.DtFunc.F(t, nil)
		}

		// fix DtOut
		if stg.Control.DtoFcn == "" {
			if stg.Control.DtOut < 1e-14 {
				stg.Control.DtOut = stg.Control.Dt
				stg.Control.DtoFunc = stg.Control.DtFunc
			} else {
				if stg.Control.DtOut < stg.Control.Dt {
					stg.Control.DtOut = stg.Control.Dt
				}
				stg.Control.DtoFunc = &fun.Cte{C: stg.Control.DtOut}
			}
		} else {
			stg.Control.DtoFunc, err = o.Functions.Get(stg.Control.DtoFcn)
			if err != nil {
				chk.Panic("%v", err)
			}
			stg.Control.DtOut = stg.Control.DtoFunc.F(t, nil)
		}

		// first stage
		if i == 0 {

			// gravity
			found := false
			for _, econd := range stg.EleConds {
				for j, key := range econd.Keys {
					if key == "g" {
						gfcn, err := o.Functions.Get(econd.Funcs[j])
						if err != nil {
							chk.Panic("ReadSim: cannot find function named %q corresponding to gravity constant @ stage 0\n%v", econd.Funcs[j], err)
						}
						o.Grav0 = gfcn.F(0, nil)
						found = true
						break
					}
				}
				if found {
					break
				}
			}
		}

		// update time
		t += stg.Control.Tf
	}

	// read materials database and initialise models
	o.Data.Wlevel = utl.Max(o.Data.Wlevel, o.MaxElev)
	o.MatModels, err = ReadMat(dir, o.Data.Matfile, o.Ndim, o.Data.Pstress, o.Data.Wlevel, o.Grav0)
	if err != nil {
		chk.Panic("loading materials and initialising models failed:\n%v", err)
	}

	// adjustable and random parameters
	o.adjmap = make(map[int]*fun.Prm)
	for _, mat := range o.MatModels.Materials {
		for _, prm := range mat.Prms {
			o.append_adjustable_parameter(prm)
		}
	}
	for _, fcn := range o.Functions {
		for _, prm := range fcn.Prms {
			o.append_adjustable_parameter(prm)
		}
	}
	err = o.AdjRandom.Init()
	if err != nil {
		chk.Panic("cannot initialise random variables:\n%v", err)
	}

	// connect dependent adjustable parameters
	var ok bool
	for _, prm := range o.AdjDependent {
		prm.Other, ok = o.adjmap[prm.Dep]
		if !ok {
			chk.Panic("cannot find dependency dep=%d of adjustable parameter", prm.Dep)
		}
	}

	// derived: liquid model
	if o.Data.LiqMat != "" {
		lmat := o.MatModels.Get(o.Data.LiqMat)
		if lmat == nil {
			chk.Panic("cannot find liquid material named %q", o.Data.LiqMat)
		}
		if lmat.Liq == nil {
			chk.Panic("liquid model in material (%q) is not available", o.Data.LiqMat)
		}
		o.LiqMdl = lmat.Liq
	}

	// derived: gas model
	if o.Data.GasMat != "" {
		gmat := o.MatModels.Get(o.Data.GasMat)
		if gmat == nil {
			chk.Panic("cannot find gas material named %g", o.Data.GasMat)
		}
		if gmat.Gas == nil {
			chk.Panic("gas model in material (%q) is not available", o.Data.GasMat)
		}
		o.GasMdl = gmat.Gas
	}

	// fix function coefficients
	for _, fcn := range o.Functions {
		for _, prm := range fcn.Prms {
			if res, found := io.Keycode(prm.Extra, "fix"); found {
				switch res {
				case "plbot":
					if o.LiqMdl == nil {
						chk.Panic("cannot fix plbot (liquid pressure at bottom of column) in %q because liquid model is not available", fcn.Name)
					}
					prm.V, _ = o.LiqMdl.Calc(0)
				case "pgbot":
					if o.GasMdl == nil {
						chk.Panic("cannot fix pgbot (gas pressure at bottom of column) in %q because gas model is not available", fcn.Name)
					}
					prm.V, _ = o.GasMdl.Calc(0)
				}
			}
		}
	}

	// results
	return &o
}

// auxiliary ///////////////////////////////////////////////////////////////////////////////////////

// Etag2data returns the ElemData corresponding to element tag
//  Note: returns nil if not found
func (o *Region) Etag2data(etag int) *ElemData {
	for _, edat := range o.ElemsData {
		if edat.Tag == etag {
			return edat
		}
	}
	return nil
}

// GetInfo returns formatted information
func (o *Simulation) GetInfo(w goio.Writer) (err error) {
	b, err := json.MarshalIndent(o, "", "  ")
	if err != nil {
		return err
	}
	_, err = w.Write(b)
	return
}

// GetEleCond returns element condition structure by giving an elem tag
//  Note: returns nil if not found
func (o Stage) GetEleCond(elemtag int) *EleCond {
	for _, ec := range o.EleConds {
		if elemtag == ec.Tag {
			return ec
		}
	}
	return nil
}

// GetNodeBc returns node boundary condition structure by giving a node tag
//  Note: returns nil if not found
func (o Stage) GetNodeBc(nodetag int) *NodeBc {
	for _, nbc := range o.NodeBcs {
		if nodetag == nbc.Tag {
			return nbc
		}
	}
	return nil
}

// GetFaceBc returns face boundary condition structure by giving a face tag
//  Note: returns nil if not found
func (o Stage) GetFaceBc(facetag int) *FaceBc {
	for _, fbc := range o.FaceBcs {
		if facetag == fbc.Tag {
			return fbc
		}
	}
	return nil
}

// extra settings //////////////////////////////////////////////////////////////////////////////////

// SetDefault sets defaults values
func (o *LinSolData) SetDefault() {
	o.Name = "umfpack"
	o.Ordering = "amf"
	o.Scaling = "rcit"
}

// SetDefault set defaults values
func (o *SolverData) SetDefault() {

	// nonlinear solver
	o.Type = "imp"
	o.NmaxIt = 20
	o.Atol = 1e-6
	o.Rtol = 1e-6
	o.FbTol = 1e-8
	o.FbMin = 1e-14
	o.NdvgMax = 20

	// Richardson's extrapolation
	o.REnssmax = 10000
	o.REatol = 1e-6
	o.RErtol = 1e-6
	o.REmfac = 0.9
	o.REmmin = 0.1
	o.REmmax = 2.0

	// transient analyses
	o.DtMin = 1e-8
	o.Theta = 0.5

	// dynamics
	o.Theta1 = 0.5
	o.Theta2 = 0.5
	o.HHTalp = 0.5

	// constants
	o.Eps = 1e-16
}

// PostProcess performs a post-processing of the just read json file
func (o *SolverData) PostProcess() {

	// coefficients for transient analyses
	if o.ThGalerkin {
		o.Theta = 2.0 / 3.0
	}
	if o.ThLiniger {
		o.Theta = 0.878
	}
	if o.ThCombo1 {
		o.Theta = 2.0 / 3.0
		o.Theta1 = 5.0 / 6.0
		o.Theta2 = 8.0 / 9.0
	}

	// iterations tolerance
	o.Itol = utl.Max(10.0*o.Eps/o.Rtol, utl.Min(0.01, math.Sqrt(o.Rtol)))
}

// adjustable parameters ///////////////////////////////////////////////////////////////////////////

// PrmAdjust adjusts parameter (random variable or not)
func (o *Simulation) PrmAdjust(adj int, val float64) {
	if prm, ok := o.adjmap[adj]; ok {
		prm.Set(val)
		return
	}
	chk.Panic("cannot adjust parameter %q", adj)
}

// PrmGetAdj gets adjustable parameter (random variable or not)
func (o *Simulation) PrmGetAdj(adj int) (val float64) {
	if prm, ok := o.adjmap[adj]; ok {
		return prm.V
	}
	chk.Panic("cannot get adjustable parameter %q", adj)
	return
}

// append_adjustable_parameter add prm to lists
func (o *Simulation) append_adjustable_parameter(prm *fun.Prm) {

	// adjustable parameter
	if prm.Adj > 0 {
		o.Adjustable = append(o.Adjustable, prm)
		o.adjmap[prm.Adj] = prm
		if prm.D != "" { // with probability distribution => random variable
			distr := rnd.GetDistribution(prm.D)
			o.AdjRandom = append(o.AdjRandom, &rnd.VarData{
				D: distr, M: prm.V, S: prm.S, Min: prm.Min, Max: prm.Max, Prm: prm,
				Key: io.Sf("adj%d", prm.Adj),
			})
		}
	}

	// adjustable parameter that depend on other adjustable parameters
	if prm.Dep > 0 {
		o.AdjDependent = append(o.AdjDependent, prm)
	}
}
