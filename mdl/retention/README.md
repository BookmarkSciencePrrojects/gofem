# Package retention implements models for liquid retention curves

## Structures and Interfaces

### Model

*Model* implements a liquid retention model (LRM)

Derivs computes (see [1] page 618):
  L  = ∂Cc/∂pc
  Lx = ∂²Cc/∂pc²
  J  = ∂Cc/∂sl
  Jx == ∂²Cc/(∂pc ∂sl)
  Jy == ∂²Cc/∂sl²

*References*:

[1] Pedroso DM (2015) A consistent u-p formulation for porous media with hysteresis.
   Int Journal for Numerical Methods in Engineering, 101(8) 606-634
   http://dx.doi.org/10.1002/nme.4808


## Models

Currently implemented models:


### BrooksCorey implements Books and Corey' model


### Lin implements a linear retetion model: sl(pc) := 1 - λ*pc


### RefM1 implements a nonlinear liquid retention model based on the concept of references [1,2,3]

*References*:

[1] Pedroso DM, Sheng D and Zhao, J (2009) The concept of reference curves for constitutive
    modelling in soil mechanics, Computers and Geotechnics, 36(1-2), 149-165,
    http://dx.doi.org/10.1016/j.compgeo.2008.01.009
[2] Pedroso DM and Williams DJ (2010) A novel approach for modelling soil-water
    characteristic curves with hysteresis, Computers and Geotechnics, 37(3), 374-380,
    http://dx.doi.org/10.1016/j.compgeo.2009.12.004
[3] Pedroso DM and Williams DJ (2011) Automatic Calibration of soil-water characteristic
    curves using genetic algorithms. Computers and Geotechnics, 38(3), 330-340,
    http://dx.doi.org/10.1016/j.compgeo.2010.12.004

### VanGen implements van Genuchten's model
