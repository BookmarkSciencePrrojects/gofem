// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/fun/dbf"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

// PlotFdata holds information to plot functions
type PlotFdata struct {
	Ti      float64  `json:"ti"`      // initial time
	Tf      float64  `json:"tf"`      // final time
	Np      int      `json:"np"`      // number of points
	Skip    []string `json:"skip"`    // skip functions
	WithTxt bool     `json:"withtxt"` // show text corresponding to initial and final points
}

// FuncData holds function definition
type FuncData struct {
	Name     string     `json:"name"`     // name of function. ex: zero, load, myfunction1, etc.
	Type     string     `json:"type"`     // type of function. ex: cte, rmp
	Prms     dbf.Params `json:"prms"`     // parameters
	PltExtra string     `json:"pltextra"` // extra arguments for plotting

	// extra data for plotting
	LabelT, LabelF, LabelG, LabelH, ArgsF, ArgsG, ArgsH string
}

// Funcs holds functions
type FuncsData []*FuncData

// Get returns function by name
func (o FuncsData) Get(name string) (fcn fun.TimeSpace, err error) {
	if name == "zero" || name == "none" {
		fcn = &fun.Zero
		return
	}
	for _, f := range o {
		if f.Name == name {
			fcn, err = fun.New(f.Type, f.Prms)
			if err != nil {
				err = chk.Err("cannot get function named %q because of the following error:\n%v", name, err)
			}
			return
		}
	}
	err = chk.Err("cannot find function named %q\n", name)
	return
}

// PlotAll plot all functions
func (o FuncsData) PlotAll(pd *PlotFdata, dirout, fnkey string) {
	for _, f := range o {
		if utl.StrIndexSmall(pd.Skip, f.Name) >= 0 {
			continue
		}
		ff, err := o.Get(f.Name)
		if err != nil {
			chk.Panic("%v", err)
		}
		if ff != nil {
			plt.Reset(false, nil)
			if pd.WithTxt {
				x := pd.Ti
				y := ff.F(x, nil)
				plt.Text(x, y, io.Sf("%g", y), nil)
				x = pd.Tf
				y = ff.F(x, nil)
				plt.Text(x, y, io.Sf("%g", y), nil)
			}
			/* TODO
			fun.PlotT(ff, "", "", pd.Ti, pd.Tf, nil, pd.Np,
				f.LabelT, f.LabelF, f.LabelG, f.LabelH, f.ArgsF, f.ArgsG, f.ArgsH)
			*/
			plt.Save(dirout, io.Sf("functions-%s-%s", fnkey, f.Name))
		}
	}
}

// auxiliary //////////////////////////////////////////////////////////////////////////////////////////

// String prints one function
func (o FuncData) String() string {
	fun.G_extraindent = "        "
	return io.Sf("    {\n      \"name\":%q, \"type\":%q, \"prms\" : [\n%v\n      ]\n    }", o.Name, o.Type, o.Prms)
}

// String prints functions
func (o FuncsData) String() string {
	if len(o) == 0 {
		return "  \"functions\" : []"
	}
	l := "  \"functions\" : [\n"
	for i, f := range o {
		if i > 0 {
			l += ",\n"
		}
		l += io.Sf("%v", f)
	}
	l += "\n  ]"
	return l
}
