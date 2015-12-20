# Gofem &ndash; Go Finite Element Method

Gofem is an implementation of the Finite Element Method in Go programming
language (golang).

## License

Unless otherwise noted, the Gofem source files are distributed
under the BSD-style license found in the LICENSE file.

## Examples

See examples here: https://github.com/cpmech/gofem/blob/master/examples/README.md

## Installation and documentation

http://cpmech.github.io/gofem

## Acknowledgements
Funding from the Australian Research Council is gratefully acknowledged.

## Subpackages are
1. ana -- analytical solutions
2. fem -- finite element method implementation; elements, solver, ...
3. inp -- input data structures; .sim (simulation), .msh (mesh), .mat (materials)
4. mdl -- material models: cnd (conductivity), fld (fluid), lrm (liquid retention), por (porous medium), sld (solids)
5. out -- output routines
5. shp -- shape (interpolation) functions; geometric elements

Additionally, 'tools' contains some auxiliary tools.
