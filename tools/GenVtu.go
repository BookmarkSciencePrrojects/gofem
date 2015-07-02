// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"bytes"
	"flag"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gofem/shp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
)

const IP_TAG_INI = 70000

var (
	ndim  int         // space dimension
	verts []*inp.Vert // all vertices
	cells []*inp.Cell // all cells
	nodes []*fem.Node // active/allocated nodes
	elems []fem.Elem  // active/allocated elements

	ipvals []map[string]float64 // [allNip][nkeys] integration points values
	exvals []map[string]float64 // [nverts][nkeys] extrapolated values
	excnts []map[string]float64 // [nverts][nkeys] extrapolation: counter

	dirout string // directory for output
	fnkey  string // filename key
	steady bool   // steady simulation

	ukeys   = []string{"ux", "uy", "uz"}                      // displacement keys
	skeys   = []string{"sx", "sy", "sz", "sxy", "syz", "szx"} // stress keys
	nwlkeys = []string{"nwlx", "nwly", "nwlz"}                // nl・wl == filter velocity keys
	plkeys  = []string{"pl"}                                  // liquid pressure keys
	pgkeys  = []string{"pg"}                                  // gas pressure keys
	flkeys  = []string{"fl"}                                  // constraint/flux/seepage face key

	is_sig     map[string]bool     // is sigma key? "sx" => true
	is_nwl     map[string]bool     // is nwl key? "nwlx" => true
	label2keys map[string][]string // maps, e.g., "u" => ukeys
)

func init() {
	is_sig = map[string]bool{"sx": true, "sy": true, "sz": true, "sxy": true, "syz": true, "szx": true}
	is_nwl = map[string]bool{"nwlx": true, "nwly": true, "nwlz": true}
	label2keys = map[string][]string{
		"u": ukeys, "sig": skeys, "nwl": nwlkeys, "ex_nwl": nwlkeys,
		"pl": plkeys, "pg": pgkeys, "fl": flkeys,
	}
}

func main() {

	// finalise analysis process and catch errors
	defer out.End()

	// input data
	simfn := "data/twoqua4.sim"
	exnwl := false
	stgidx := 0

	// parse flags
	flag.Parse()
	if len(flag.Args()) > 0 {
		simfn = flag.Arg(0)
	}
	if len(flag.Args()) > 1 {
		exnwl = io.Atob(flag.Arg(1))
	}
	if len(flag.Args()) > 2 {
		stgidx = io.Atoi(flag.Arg(2))
	}

	// check extension
	if io.FnExt(simfn) == "" {
		simfn += ".sim"
	}

	// print input data
	io.Pf("\nInput data\n")
	io.Pf("==========\n")
	io.Pf("  simfn   = %30s // simulation filename\n", simfn)
	io.Pf("  exnwl   = %30v // extrapolate nwl\n", exnwl)
	io.Pf("  stgidx  = %30v // stage index\n", stgidx)
	io.Pf("\n")

	// start analysis process
	out.Start(simfn, stgidx, 0)

	// global variables
	ndim = out.Dom.Msh.Ndim
	verts = out.Dom.Msh.Verts
	cells = out.Dom.Msh.Cells
	nodes = out.Dom.Nodes
	elems = out.Dom.Elems
	dirout = fem.Global.Sim.Data.DirOut
	fnkey = fem.Global.Sim.Data.FnameKey
	steady = fem.Global.Sim.Data.Steady

	// flags
	has_u := out.Dom.YandC["ux"]
	has_pl := out.Dom.YandC["pl"]
	has_pg := out.Dom.YandC["pg"]
	has_sig := out.Ipkeys["sx"]
	has_nwl := out.Ipkeys["nwlx"]
	has_p := has_pl || has_pg
	lbb := has_u && has_p
	if fem.Global.Sim.Data.NoLBB {
		lbb = false
	}

	// buffers
	pvd := make(map[string]*bytes.Buffer)
	geo := make(map[string]*bytes.Buffer)
	vtu := make(map[string]*bytes.Buffer)
	if _, ok := out.Dom.YandC["ux"]; ok {
		pvd["u"] = new(bytes.Buffer)
		geo["u"] = new(bytes.Buffer)
		vtu["u"] = new(bytes.Buffer)
	}
	if _, ok := out.Dom.YandC["fl"]; ok {
		pvd["fl"] = new(bytes.Buffer)
		geo["fl"] = new(bytes.Buffer)
		vtu["fl"] = new(bytes.Buffer)
	}
	if _, ok := out.Dom.YandC["pl"]; ok {
		pvd["pl"] = new(bytes.Buffer)
		geo["pl"] = new(bytes.Buffer)
		vtu["pl"] = new(bytes.Buffer)
	}
	if _, ok := out.Dom.YandC["pg"]; ok {
		pvd["pg"] = new(bytes.Buffer)
		geo["pg"] = new(bytes.Buffer)
		vtu["pg"] = new(bytes.Buffer)
	}
	if len(out.Ipkeys) > 0 {
		pvd["ips"] = new(bytes.Buffer)
		geo["ips"] = new(bytes.Buffer)
		vtu["ips"] = new(bytes.Buffer)
	}
	if exnwl {
		pvd["ex_nwl"] = new(bytes.Buffer)
		geo["ex_nwl"] = new(bytes.Buffer)
		vtu["ex_nwl"] = new(bytes.Buffer)
	}

	// headers
	for _, b := range pvd {
		pvd_header(b)
	}

	// process results
	for tidx, t := range out.Sum.OutTimes {

		// input results into domain
		if !out.Dom.In(out.Sum, tidx, true) {
			chk.Panic("cannot load results into domain; please check log file")
		}

		// message
		io.PfWhite("time     = %g\r", t)

		// generate topology
		if tidx == 0 {
			for label, b := range geo {
				topology(b, label == "ips", lbb)
			}

			// allocate integration points values
			ipvals = make([]map[string]float64, len(out.Ipoints))
			for i, _ := range out.Ipoints {
				ipvals[i] = make(map[string]float64)
			}
		}

		// get integration points values @ time t
		for i, p := range out.Ipoints {
			vals := p.Calc(out.Dom.Sol)
			for key, val := range vals {
				ipvals[i][key] = val
			}
		}

		// compute extrapolated values
		if exnwl {
			compute_extrapolated_values()
		}

		// for each data buffer
		for label, b := range vtu {

			// reset buffer
			b.Reset()

			// points data
			if label == "ips" {
				pdata_open(b)
				if has_sig {
					pdata_write(b, "sig", skeys, true)
				}
				if has_nwl {
					pdata_write(b, "nwl", nwlkeys, true)
				}
				for key, _ := range out.Ipkeys {
					if !is_sig[key] && !is_nwl[key] {
						pdata_write(b, key, []string{key}, true)
					}
				}
				pdata_close(b)
			} else {
				pdata_open(b)
				pdata_write(b, label, label2keys[label], false)
				pdata_close(b)
			}

			// cells data
			cdata_write(b, label == "ips")

			// write vtu file
			vtu_write(geo[label], b, tidx, label)
		}

		// pvd
		for label, b := range pvd {
			pvd_line(b, tidx, t, label)
		}
	}

	// write pvd files
	for label, b := range pvd {
		pvd_write(b, label)
	}
}

// headers and footers ///////////////////////////////////////////////////////////////////////////////

func pvd_header(buf *bytes.Buffer) {
	if buf == nil {
		return
	}
	io.Ff(buf, "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n<Collection>\n")
}

func pvd_line(buf *bytes.Buffer, tidx int, time float64, label string) {
	if buf == nil {
		return
	}
	io.Ff(buf, "<DataSet timestep=\"%23.15e\" file=\"%s_%06d_%s.vtu\" />\n", time, fnkey, tidx, label)
}

func pvd_write(buf *bytes.Buffer, label string) {
	if buf == nil {
		return
	}
	io.Ff(buf, "</Collection>\n</VTKFile>")
	io.WriteFileV(io.Sf("%s/%s_%s.pvd", dirout, fnkey, label), buf)
}

func vtu_write(geo, dat *bytes.Buffer, tidx int, label string) {
	if geo == nil || dat == nil {
		return
	}
	nv := len(verts)
	nc := len(elems)
	if label == "ips" {
		nv = len(out.Ipoints)
		nc = nv
	}
	var hdr, foo bytes.Buffer
	io.Ff(&hdr, "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n")
	io.Ff(&hdr, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nv, nc)
	io.Ff(&foo, "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n")
	io.WriteFile(io.Sf("%s/%s_%06d_%s.vtu", dirout, fnkey, tidx, label), &hdr, geo, dat, &foo)
}

// topology ////////////////////////////////////////////////////////////////////////////////////////

func topology(buf *bytes.Buffer, ips, lbb bool) {
	if buf == nil {
		return
	}

	// coordinates
	io.Ff(buf, "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n")
	var z float64
	if ips {
		for _, p := range out.Ipoints {
			if ndim == 3 {
				z = p.X[2]
			}
			io.Ff(buf, "%23.15e %23.15e %23.15e ", p.X[0], p.X[1], z)
		}
	} else {
		for _, v := range verts {
			if ndim == 3 {
				z = v.C[2]
			}
			io.Ff(buf, "%23.15e %23.15e %23.15e ", v.C[0], v.C[1], z)
		}
	}
	io.Ff(buf, "\n</DataArray>\n</Points>\n")

	// connectivities
	io.Ff(buf, "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")
	if ips {
		for i, _ := range out.Ipoints {
			io.Ff(buf, "%d ", i)
		}
	} else {
		for _, e := range elems {
			cell := cells[e.Id()]
			_, nverts := get_cell_info(cell.Type, lbb)
			if cell.Type == "joint" {
				nverts = len(cell.Verts)
			}
			for j := 0; j < nverts; j++ {
				io.Ff(buf, "%d ", cell.Verts[j])
			}
		}
	}

	// offsets of active elements
	io.Ff(buf, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")
	var offset int
	if ips {
		for range out.Ipoints {
			offset += 1
			io.Ff(buf, "%d ", offset)
		}
	} else {
		for _, e := range elems {
			cell := cells[e.Id()]
			_, nverts := get_cell_info(cell.Type, lbb)
			if cell.Type == "joint" {
				nverts = len(cell.Verts)
			}
			offset += nverts
			io.Ff(buf, "%d ", offset)
		}
	}

	// types of active elements
	io.Ff(buf, "\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n")
	if ips {
		for range out.Ipoints {
			io.Ff(buf, "%d ", shp.VTK_VERTEX)
		}
	} else {
		for _, e := range elems {
			cell := cells[e.Id()]
			ctype, _ := get_cell_info(cell.Type, lbb)
			vtk := shp.GetVtkCode(ctype)
			if ctype == "joint" {
				vtk = shp.VTK_POLY_VERTEX
			}
			if vtk < 0 {
				chk.Panic("cannot handle cell type %q", cell.Type)
			}
			io.Ff(buf, "%d ", vtk)
		}
	}
	io.Ff(buf, "\n</DataArray>\n</Cells>\n")
	return
}

// points data /////////////////////////////////////////////////////////////////////////////////////

func pdata_open(buf *bytes.Buffer) {
	if buf == nil {
		return
	}
	io.Ff(buf, "<PointData Scalars=\"TheScalars\">\n")
}

func pdata_close(buf *bytes.Buffer) {
	if buf == nil {
		return
	}
	io.Ff(buf, "</PointData>\n")
}

func get_zeros(n int) (l string) {
	for i := 0; i < n; i++ {
		l += "0 "
	}
	return
}

func pdata_write(buf *bytes.Buffer, label string, keys []string, ips bool) {
	if buf == nil {
		return
	}

	nkeys := len(keys)
	zeros := get_zeros(nkeys)
	extrap := false
	if len(label) > 2 {
		extrap = label[:2] == "ex"
	}

	io.Ff(buf, "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\">\n", label, nkeys)
	if ips {
		// loop over integration points
		for i, _ := range out.Ipoints {
			l := ""
			for _, key := range keys {
				if val, ok := ipvals[i][key]; ok {
					l += io.Sf("%23.15e ", val)
				} else {
					l += "0 "
				}
			}
			if l == "" {
				l = zeros
			}
			io.Ff(buf, l)
		}
	} else if extrap {
		// loop over vertices => with extrapolated values
		for _, v := range verts {
			n := out.Dom.Vid2node[v.Id]
			l := zeros
			if n != nil {
				l = ""
				for _, key := range keys {
					if val, ok := exvals[v.Id][key]; ok {
						l += io.Sf("%23.15e ", val/excnts[v.Id][key])
					} else {
						l += "0 "
					}
				}
			}
			io.Ff(buf, l)
		}
	} else {
		// loop over vertices
		Y := out.Dom.Sol.Y
		switch label {
		case "v":
			Y = out.Dom.Sol.Dydt
		case "a":
			Y = out.Dom.Sol.D2ydt2
		}
		for _, v := range verts {
			n := out.Dom.Vid2node[v.Id]
			l := zeros
			if n != nil {
				l = ""
				for _, key := range keys {
					eq := n.GetEq(key)
					if eq >= 0 {
						l += io.Sf("%23.15e ", Y[eq])
					} else {
						l += "0 "
					}
				}
			}
			io.Ff(buf, l)
		}
	}
	io.Ff(buf, "\n</DataArray>\n")
}

func cdata_write(buf *bytes.Buffer, ips bool) {
	if buf == nil {
		return
	}

	// open
	io.Ff(buf, "<CellData Scalars=\"TheScalars\">\n")

	// ids
	io.Ff(buf, "<DataArray type=\"Float64\" Name=\"eid\" NumberOfComponents=\"1\" format=\"ascii\">\n")
	if ips {
		for _, p := range out.Ipoints {
			io.Ff(buf, "%d ", p.Eid)
		}
	} else {
		for _, e := range elems {
			io.Ff(buf, "%d ", e.Id())
		}
	}

	// cells positive tags
	io.Ff(buf, "\n</DataArray>\n<DataArray type=\"Float64\" Name=\"tag\" NumberOfComponents=\"1\" format=\"ascii\">\n")
	if ips {
		for _, p := range out.Ipoints {
			ptag := IP_TAG_INI + iabs(cells[p.Eid].Tag)
			io.Ff(buf, "%d ", ptag)
		}
	} else {
		for _, e := range elems {
			ptag := iabs(cells[e.Id()].Tag)
			io.Ff(buf, "%d ", ptag)
		}
	}

	// close
	io.Ff(buf, "\n</DataArray>\n</CellData>\n")
}

func iabs(val int) int {
	if val < 0 {
		return -val
	}
	return val
}

func get_cell_info(ctype string, lbb bool) (ctypeNew string, nverts int) {
	ctypeNew = ctype
	if ctypeNew == "qua9" {
		ctypeNew = "qua8"
	}
	if lbb {
		ctypeNew = shp.GetBasicType(ctypeNew)
	}
	nverts = shp.GetNverts(ctypeNew)
	return
}

func compute_extrapolated_values() {

	// allocate structures for extrapolation
	if len(exvals) == 0 {
		exvals = make([]map[string]float64, len(verts))
		excnts = make([]map[string]float64, len(verts))
	}

	// clear previous values
	for _, v := range verts {
		exvals[v.Id] = make(map[string]float64)
		excnts[v.Id] = make(map[string]float64)
	}

	// keys
	keys := []string{"nwlx", "nwly"}
	if ndim == 2 {
		keys = []string{"nwlx", "nwly", "nwlz"}
	}

	// loop over elements
	for _, ele := range elems {

		// get shape and integration points from known elements
		var sha *shp.Shape
		var ips []*shp.Ipoint
		switch e := ele.(type) {
		case *fem.ElemP:
			sha = e.Shp
			ips = e.IpsElem
		case *fem.ElemU:
			sha = e.Shp
			ips = e.IpsElem
		case *fem.ElemUP:
			sha = e.U.Shp
			ips = e.U.IpsElem
		}
		if sha == nil {
			chk.Panic("cannot get shape structure from element")
		}

		// compute Extrapolator matrix
		Emat := la.MatAlloc(sha.Nverts, len(ips))
		err := sha.Extrapolator(Emat, ips)
		if err != nil {
			chk.Panic("cannot compute extrapolator matrix: %v", err)
		}

		// perform extrapolation
		cell := cells[ele.Id()]
		ids := out.Cid2ips[ele.Id()]
		for i := 0; i < sha.Nverts; i++ {
			v := cell.Verts[i]
			for _, key := range keys {
				for j := 0; j < len(ips); j++ {
					exvals[v][key] += Emat[i][j] * ipvals[ids[j]][key]
				}
				excnts[v][key] += 1
			}
		}
	}
}
