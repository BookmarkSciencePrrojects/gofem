// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"bytes"

	"github.com/cpmech/gofem/ele"
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gofem/shp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

const IP_TAG_INI = 70000

var (
	ndim  int           // space dimension
	verts []*inp.Vert   // all vertices
	cells []*inp.Cell   // all cells
	nodes []*fem.Node   // active/allocated nodes
	elems []ele.Element // active/allocated elements

	dirout string // directory for output
	fnkey  string // filename key
	steady bool   // steady simulation

	ukeys   = []string{"ux", "uy", "uz"}                      // displacement keys
	skeys   = []string{"sx", "sy", "sz", "sxy", "syz", "szx"} // stress keys
	nwlkeys = []string{"nwlx", "nwly", "nwlz"}                // nl・wl == liquid filter velocity keys
	nwgkeys = []string{"nwgx", "nwgy", "nwgz"}                // ng・wg == gas filter velocity keys
	plkeys  = []string{"pl"}                                  // liquid pressure keys
	pgkeys  = []string{"pg"}                                  // gas pressure keys
	flkeys  = []string{"fl"}                                  // constraint/flux/seepage face key

	is_sig     map[string]bool     // is sigma key? "sx" => true
	is_nwl     map[string]bool     // is nwl key? "nwlx" => true
	is_nwg     map[string]bool     // is nwg key? "nwgx" => true
	label2keys map[string][]string // maps, e.g., "u" => ukeys
)

func init() {
	is_sig = map[string]bool{"sx": true, "sy": true, "sz": true, "sxy": true, "syz": true, "szx": true}
	is_nwl = map[string]bool{"nwlx": true, "nwly": true, "nwlz": true}
	is_nwg = map[string]bool{"nwgx": true, "nwgy": true, "nwgz": true}
	label2keys = map[string][]string{
		"u": ukeys, "sig": skeys, "nwl": nwlkeys, "ex_nwl": nwlkeys, "nwg": nwgkeys, "ex_nwg": nwgkeys,
	}
}

func main() {

	// catch errors
	defer func() {
		if err := recover(); err != nil {
			io.PfRed("ERROR: %v\n", err)
		}
	}()

	// input data
	simfn, _ := io.ArgToFilename(0, "data/twoqua4", ".sim", true)
	exnwl := io.ArgToBool(1, false)
	exnwg := io.ArgToBool(2, false)
	stgidx := io.ArgToInt(3, 0)
	v3beam := io.ArgToBool(4, false)
	io.Pf("\n%s\n", io.ArgsTable("INPUT ARGUMENTS",
		"simulation filename", "simfn", simfn,
		"extrapolate nwl", "exnwl", exnwl,
		"extrapolate nwg", "exnwg", exnwg,
		"stage index", "stgidx", stgidx,
		"show v3 of beams", "v3beam", v3beam,
	))

	// start analysis process
	out.Start(simfn, stgidx, 0)

	// global variables
	ndim = out.Dom.Msh.Ndim
	verts = out.Dom.Msh.Verts
	cells = out.Dom.Msh.Cells
	nodes = out.Dom.Nodes
	elems = out.Dom.Elems
	dirout = out.Dom.Sim.DirOut
	fnkey = out.Dom.Sim.Key
	steady = out.Dom.Sim.Data.Steady

	// flags
	has_u := out.Dom.YandC["ux"]
	has_pl := out.Dom.YandC["pl"]
	has_pg := out.Dom.YandC["pg"]
	has_sig := out.Ipkeys["sx"]
	has_nwl := out.Ipkeys["nwlx"]
	has_nwg := out.Ipkeys["nwgx"]
	has_p := has_pl || has_pg
	lbb := has_u && has_p
	if out.Dom.Sim.Data.NoLBB {
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
	for ykey, _ := range out.Dom.Dof2Tnum {
		if ykey == "ux" || ykey == "uy" || ykey == "uz" {
			continue
		}
		pvd[ykey] = new(bytes.Buffer)
		geo[ykey] = new(bytes.Buffer)
		vtu[ykey] = new(bytes.Buffer)
		label2keys[ykey] = []string{ykey}
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
	if exnwg {
		pvd["ex_nwg"] = new(bytes.Buffer)
		geo["ex_nwg"] = new(bytes.Buffer)
		vtu["ex_nwg"] = new(bytes.Buffer)
	}

	// extrapolated values keys
	var extrap_keys []string
	if exnwl && has_nwl {
		extrap_keys = append(extrap_keys, nwlkeys[:ndim]...)
	}
	if exnwg && has_nwg {
		extrap_keys = append(extrap_keys, nwgkeys[:ndim]...)
	}

	// headers
	for _, b := range pvd {
		pvd_header(b)
	}

	// process results
	for tidx, t := range out.Sum.OutTimes {

		// input results into domain
		err := out.Dom.Read(out.Sum, tidx, 0, true)
		if err != nil {
			chk.Panic("cannot load results into domain\n%v", err)
		}

		// message
		io.PfWhite("time     = %g\r", t)

		// generate topology
		if tidx == 0 {
			for label, b := range geo {
				topology(b, label == "ips", lbb, v3beam)
			}
		}

		// current values @ ips
		for _, element := range out.ElemOutIps {
			ipids := out.Cid2ips[element.Id()]
			allvals := ele.NewIpsMap()
			element.OutIpVals(allvals, out.Dom.Sol)
			for key, vals := range *allvals {
				for i, ipid := range ipids {
					out.Ipoints[ipid].Vals[key] = vals[i]
				}
			}
		}

		// compute extrapolated values
		if len(extrap_keys) > 0 {
			out.ComputeExtrapolatedValues(extrap_keys)
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
				if has_nwg {
					pdata_write(b, "nwg", nwgkeys, true)
				}
				for key, _ := range out.Ipkeys {
					if !is_sig[key] && !is_nwl[key] && !is_nwg[key] {
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

func topology(buf *bytes.Buffer, ips, lbb, v3beam bool) {
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
			nverts, _ := cell.GetVtkInfo(lbb, v3beam)
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
			nverts, _ := cell.GetVtkInfo(lbb, v3beam)
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
			_, vtkcode := cell.GetVtkInfo(lbb, v3beam)
			if vtkcode < 0 {
				chk.Panic("cannot handle cell type %q", cell.Shp.Type)
			}
			io.Ff(buf, "%d ", vtkcode)
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
		for _, ip := range out.Ipoints {
			l := ""
			for _, key := range keys {
				if val, ok := ip.Vals[key]; ok {
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
					if val, ok := out.ExVals[v.Id][key]; ok {
						l += io.Sf("%23.15e ", val)
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

	// positive tags
	if !ips {
		io.Ff(buf, "<DataArray type=\"Int32\" Name=\"tag\" NumberOfComponents=\"1\" format=\"ascii\">\n")
		for _, v := range verts {
			io.Ff(buf, "%d ", iabs(v.Tag))
		}
		io.Ff(buf, "\n</DataArray>\n")
	}
}

func cdata_write(buf *bytes.Buffer, ips bool) {
	if buf == nil {
		return
	}

	// open
	io.Ff(buf, "<CellData Scalars=\"TheScalars\">\n")

	// ids
	io.Ff(buf, "<DataArray type=\"Int32\" Name=\"eid\" NumberOfComponents=\"1\" format=\"ascii\">\n")
	if ips {
		for _, p := range out.Ipoints {
			io.Ff(buf, "%d ", p.Cid)
		}
	} else {
		for _, e := range elems {
			io.Ff(buf, "%d ", e.Id())
		}
	}

	// cells positive tags
	io.Ff(buf, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"tag\" NumberOfComponents=\"1\" format=\"ascii\">\n")
	if ips {
		for _, p := range out.Ipoints {
			ptag := IP_TAG_INI + iabs(cells[p.Cid].Tag)
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
