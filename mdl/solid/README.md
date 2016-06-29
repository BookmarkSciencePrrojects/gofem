# Package sld implements models for solids based on continuum mechanics

## Structures and Interfaces

*Model* defines the interface for solid models

*Driver* run simulations with constitutive models for solids

*Plotter* assists on plotting numerical results

*State* holds all continuum mechanics data, including for updating the state



## Models

*CamClayMod* implements the modified CamClay model

*DruckerPrager* implements Drucker-Prager plasticity model

*SmallElasticity* implements linear/non-linear elasticity for small strain analyses

*HyperElast1* implements a nonlinear hyperelastic model for powders and porous media

*LinElast* implements a linear elastic model

*Ogden* implements a linear elastic model

*OnedLinElast* implements a linear elastic model for 1D elements

*RjointM1* implements a 1D plasticity model for rod-joints (links/interface)

*SmpInvs* implements a model with SMP invariants similar to Drucker-Prager model
