# Package fld implements models for fluid density

## Structures and Interfaces

### Model

*Model* implements a model to compute pressure (p) and intrinsic density (R) of a fluid
along a column with gravity (g).

The model is:
  R(p) = R0 + C・(p - p0)   thus   dR/dp = C
