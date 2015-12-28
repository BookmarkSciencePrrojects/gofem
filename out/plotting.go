// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"math"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

// PltEntity stores all data for a plot entity (X vs Y)
type PltEntity struct {
	Alias string    // alias
	X     []float64 // x-values
	Y     []float64 // y-values
	Xlbl  string    // horizontal axis label (raw; e.g. "t")
	Ylbl  string    // vertical axis label (raw; e.g. "pl")
	Style plt.Fmt   // style
}

// SplotDat stores all data for one subplot
type SplotDat struct {
	Id      string       // unique identifier
	Title   string       // title of subplot
	Topts   string       // title options
	Xscale  float64      // x-axis scale
	Yscale  float64      // y-axis scale
	Xrange  []float64    // x range
	Yrange  []float64    // x range
	Xlbl    string       // x-axis label (formatted; e.g. "$t$")
	Ylbl    string       // y-axis label (formatted; e.g. "$p_{\ell}$")
	GllArgs string       // extra arguments for Gll such as leg_out
	Data    []*PltEntity // data and styles to be plotted
}

// Splot activates a new subplot window
func Splot(id, splotTitle string) {
	s := &SplotDat{Id: id, Title: splotTitle}
	Splots = append(Splots, s)
	Csplot = s
}

// SplotConfig configures units and scales of axes
func SplotConfig(xunit, yunit string, xscale, yscale float64) {
	if Csplot != nil {
		var xlabel, ylabel string
		if len(Csplot.Data) > 0 {
			xlabel = Csplot.Data[0].Xlbl
			ylabel = Csplot.Data[0].Ylbl
		}
		Csplot.Xlbl = GetTexLabel(xlabel, xunit)
		Csplot.Ylbl = GetTexLabel(ylabel, yunit)
		Csplot.Xscale = xscale
		Csplot.Yscale = yscale
	}
}

// Plot plots data
//  xHandle -- can be a string, e.g. "t" or a slice, e.g. pc = []float64{0, 1, 2}
//  yHandle -- can be a string, e.g. "pl" or a slice, e.g. sl = []float64{0, 1, 2}
//  alias   -- alias such as "centre"
//  fm      -- formatting codes; e.g. plt.Fmt{C:"blue", L:"label"}
//  idxI    -- index of time; use -1 for all times
func Plot(xHandle, yHandle interface{}, alias string, fm plt.Fmt, idxI int) {
	var e PltEntity
	e.Alias = alias
	e.Style = fm
	e.X, e.Xlbl = get_vals_and_labels(xHandle, yHandle, alias, idxI)
	e.Y, e.Ylbl = get_vals_and_labels(yHandle, xHandle, alias, idxI)
	if len(e.X) != len(e.Y) {
		chk.Panic("lengths of x- and y-series are different. len(x)=%d, len(y)=%d, x=%v, y=%v", len(e.X), len(e.Y), xHandle, yHandle)
	}
	if Csplot == nil {
		Splot(io.Sf("%d", len(Splots)), "")
	}
	Csplot.Data = append(Csplot.Data, &e)
	SplotConfig("", "", 1, 1)
}

// Draw draws or save figure with plot
//  dirout -- directory to save figure
//  fname  -- file name; e.g. myplot.eps or myplot.png. Use "" to show figure instead
//  nr     -- number of rows. Use -1 to compute best value
//  nc     -- number of columns. Use -1 to compute best value
//  split  -- split subplots into separated figures
//  extra  -- is called just after Subplot command and before any plotting
func Draw(dirout, fname string, nr, nc int, split bool, extra func(id string)) {
	var fnk string // filename key
	var ext string // extension
	if fname != "" {
		fnk = io.FnKey(fname)
		ext = io.FnExt(fname)
	}
	nplots := len(Splots)
	if nr < 0 || nc < 0 {
		nr, nc = utl.BestSquare(nplots)
	}
	for k := 0; k < nplots; k++ {
		spl := Splots[k]
		if !split {
			plt.Subplot(nr, nc, k+1)
		}
		if extra != nil {
			extra(spl.Id)
		}
		if spl.Title != "" {
			plt.Title(spl.Title, spl.Topts)
		}
		for _, d := range spl.Data {
			if d.Style.L == "" {
				d.Style.L = d.Alias
			}
			x, y := d.X, d.Y
			if math.Abs(spl.Xscale) > 0 {
				x = make([]float64, len(d.X))
				la.VecCopy(x, spl.Xscale, d.X)
			}
			if math.Abs(spl.Yscale) > 0 {
				y = make([]float64, len(d.Y))
				la.VecCopy(y, spl.Yscale, d.Y)
			}
			plt.Plot(x, y, d.Style.GetArgs("clip_on=0"))
		}
		plt.Gll(spl.Xlbl, spl.Ylbl, spl.GllArgs)
		if len(spl.Xrange) == 2 {
			plt.AxisXrange(spl.Xrange[0], spl.Xrange[1])
		}
		if len(spl.Yrange) == 2 {
			plt.AxisYrange(spl.Yrange[0], spl.Yrange[1])
		}
		if split {
			savefig(dirout, fnk, ext, spl.Id)
			plt.Clf()
		}
	}
	if !split && fname != "" {
		savefig(dirout, fnk, ext, "")
	}
	if fname == "" {
		plt.Show()
	}
}

// auxiliary /////////////////////////////////////////////////////////////////////////////////////////

func savefig(dirout, fnk, ext, id string) {
	fn := fnk + ext
	if id != "" {
		fn = fnk + "_" + id + ext
	}
	if dirout == "" {
		plt.Save(fn)
	} else {
		plt.SaveD(dirout, fn)
	}
}

func get_vals_and_labels(handle, otherHandle interface{}, alias string, idxI int) ([]float64, string) {
	otherKey := "any"
	if key, ok := otherHandle.(string); ok {
		otherKey = key
	}
	switch hnd := handle.(type) {
	case []float64:
		return hnd, io.Sf("%s-type", alias)
	case string:
		switch hnd {
		case "t":
			return Times, "t"
		case "x":
			xcoords, _, _ := GetXYZ(otherKey, alias)
			return xcoords, "x"
		case "y":
			_, ycoords, _ := GetXYZ(otherKey, alias)
			return ycoords, "y"
		case "z":
			_, _, zcoords := GetXYZ(otherKey, alias)
			return zcoords, "z"
		case "dist":
			return GetDist(otherKey, alias), "dist"
		}
		return GetRes(hnd, alias, idxI), hnd
	}
	chk.Panic("cannot get values slice with handle = %v", handle)
	return nil, ""
}
