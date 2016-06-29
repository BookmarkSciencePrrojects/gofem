# Gofem &ndash; Go Finite Element Method

Gofem (Go Finite Element Method) is an implementation of the finite element method (FEM) in Go
language for applications in solid mechanics. The code aims to be as general as possible and has a
focus on porous media mechanics. Nonetheless, classical plasticity and the solution of multi-physics
coupled problems are also targeted by Gofem. Efficiency is a goal as long as the quality of code and
code maintenance is not penalised. The computational efficiency is achieved by parallel computing
using message passage interface (MPI). Several unit tests are employed for every detail of the code
and its usage aims to be comprehensive. Gofem depends on the Go Scientific Library (Gosl) and was
developed for obtaining the results presented in a number of journal papers, including [1-5].



## Content

1.  ana              &ndash; analytical solutions
2.  shp              &ndash; shape (interpolation) structures and quadrature points
3.  mdl/generic      &ndash; generic models (placeholder for parameters set)
4.  mdl/solid        &ndash; models for solids
5.  mdl/fluid        &ndash; models for fluids (liquid / gas)
6.  mdl/conduct      &ndash; models for liquid conductivity within porous media
7.  mdl/retention    &ndash; models for liquid retention within porous media
8.  mdl/diffusion    &ndash; models for diffusion applications
9.  mdl/thermomech   &ndash; models for thermo-mechanical applications
10. mdl/porous       &ndash; models for porous media (TPM-based)
11. inp              &ndash; input data (.sim = simulation, .mat = materials, .msh = meshes)
12. ele              &ndash; finite elements
13. ele/solid        &ndash; elements for solid mechanics
14. ele/seepage      &ndash; elements for seepage problems (with liquid and/or gases)
15. ele/diffusion    &ndash; elements for diffusion(-like) problems
16. ele/thermomech   &ndash; elements for thermo-mechanical applications
17. ele/porous       &ndash; elements for porous media simulations (TPM)
18. fem              &ndash; finite element method (main, domain, solver, ...)
19. tests            &ndash; (unit) tests of complete simulations
20. tests/solid      &ndash; tests of solid mechanics applications
21. tests/seepage    &ndash; tests of seepage problems
22. tests/diffusion  &ndash; tests of diffusion problems
23. tests/thermomech &ndash; tests of thermo-mechanical applications
24. tests/porous     &ndash; tests of porous media simulations
25. out              &ndash; output routines (post-processing and plotting)



## Examples

See examples here: https://github.com/cpmech/gofem/blob/master/examples/README.md



## Installation and documentation

```
mkdir -p $GOPATH/src/github.com/cpmech
cd $GOPATH/src/github.com/cpmech
git clone https://github.com/cpmech/gofem.git
cd gofem
./all.bash
```

See http://cpmech.github.io/gofem for more details.



## Acknowledgements
Funding from the Australian Research Council is gratefully acknowledged.

Additionally, 'tools' contains some auxiliary tools.



## References

1. Pedroso DM (2015) A consistent u-p formulation for porous media with hysteresis. Int Journal for Numerical Methods in Engineering, 101(8) 606-634 http://dx.doi.org/10.1002/nme.4808
2. Pedroso DM (2015) A solution to transient seepage in unsaturated porous media. Computer Methods in Applied Mechanics and Engineering, 285 791-816 http://dx.doi.org/10.1016/j.cma.2014.12.009
3. Pedroso DM, Sheng D and Zhao, J (2009) The concept of reference curves for constitutive modelling in soil mechanics, Computers and Geotechnics, 36(1-2), 149-165, http://dx.doi.org/10.1016/j.compgeo.2008.01.009
4. Pedroso DM and Williams DJ (2010) A novel approach for modelling soil-water characteristic curves with hysteresis, Computers and Geotechnics, 37(3), 374-380, http://dx.doi.org/10.1016/j.compgeo.2009.12.004
5. Pedroso DM and Williams DJ (2011) Automatic Calibration of soil-water characteristic curves using genetic algorithms. Computers and Geotechnics, 38(3), 330-340, http://dx.doi.org/10.1016/j.compgeo.2010.12.004



## License

Unless otherwise noted, the Gofem source files are distributed
under the BSD-style license found in the LICENSE file.
