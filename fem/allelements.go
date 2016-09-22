// Copyright 2016 The Gofem Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/ele/diffusion"
	"github.com/cpmech/gofem/ele/porous"
	"github.com/cpmech/gofem/ele/seepage"
	"github.com/cpmech/gofem/ele/solid"
	"github.com/cpmech/gofem/ele/thermomech"
)

// enforce loading of all elements
func init() {
	_ = diffusion.Diffusion{}
	_ = porous.SolidLiquid{}
	_ = seepage.Liquid{}
	_ = solid.Solid{}
	_ = thermomech.SolidThermal{}
}
