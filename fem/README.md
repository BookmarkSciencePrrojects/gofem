# Package fem implements the FEM solver

## Structures and Interfaces

1. *Main* holds all data for a simulation using the finite element method
2. *Solver* Solver implements the actual solver (time loop)
3. *Summary* records summary of outputs
4. *EssentialBc* holds information about essential bounday conditions such as constrained nodes
5. *PtNaturalBc* holds information on point natural boundary conditions such as prescribed forces or fluxes) at nodes

## Solvers

1. *SolverLinearImplicit* solves **linear** FEM problem using an implicit procedure
2. *Implicit* solves FEM problem using an implicit procedure (with Newthon-Raphson method)
3. *RichardsonExtrap* solves FEM problem implicitely and with Richardson's extrapolation
