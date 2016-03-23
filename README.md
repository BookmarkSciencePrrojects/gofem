# Gofem &ndash; Go Finite Element Method

Gofem is an implementation of the Finite Element Method in Go programming
language (golang).

## License

Unless otherwise noted, the Gofem source files are distributed
under the BSD-style license found in the LICENSE file.

## Examples

See examples here: https://github.com/cpmech/gofem/blob/master/examples/README.md

## News

2016/03/23 -- v3.0 was tagged as stable (new stable-v3.0 branch created)

## Installation and documentation (for developers)

http://cpmech.github.io/gofem

## Acknowledgements
Funding from the Australian Research Council is gratefully acknowledged.

## Subpackages
1.  ana     -- analytical solutions
2.  shp     -- shape (interpolation) structures and quadrature points
3.  inp     -- input data structures. simulation, materials, meshes
4.  mdl/sld -- models for solids
5.  mdl/fld -- models for fluids
6.  mdl/cnd -- models for liquid/gas conductivity in porous media
7.  mdl/lrm -- models for liquid retention in porous media
8.  mdl/por -- models for porous media
9.  fem     -- finite element method (elements, solver, ...)
10. out     -- output: analyses of results and plotting

Additionally, 'tools' contains some auxiliary tools.
