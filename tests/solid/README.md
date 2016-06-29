# Tests: Solid

Tests with elements for solid mechanics problems.

## Beam Element (frames)

1. beam01a. check DOFs
2. beam01b. simply supported
3. beam02. cantilever
4. beam03. small frame
5. beam04. 3D beam (bh414)
6. beam04. 3D frame

## Bhatti's Book

*Reference*

Bhatti MA (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.

1. bh16a. bracket. check DOFs
2. bh16b. bracket. run
3. bh14a. using RunAll
4. bh14b. truss. using SolveOneStage
5. bh14c. truss. call to PrmAdjust
6. bh14d. truss. using go-routines
7. bh14erod. truss. using ElasticRod

## Beam-Joint (compression) Element

1. bjoint01a. beam joint compatible. static. check
2. bjoint01b. beam joint compatible. static. run
3. bjoint02. beam joint compatible. pull-out. run
4. bjoint03. beam joint compatible. static. run
5. bjoint04. beam joint compatible. displace beam. run

## NURBS. Isogeometric analysis

1. nurbs01. square with initial stress
2. nurbs02. square with initial stress. run
3. nurbs03. ini stress free square

## Rod-Joint Element
1. rjoint01. curved line in 3D

## Rod Element (trusses)

1. bridge01a. simple bridge section
2. bridge01. simple bridge section. ElastRod

## Smith, Griffiths and Margetts' Book

*Reference*

Smith IM, Griffiths DV and Margetts L (2014) Programming the Finite Element Method, 5th Edition, Wiley, 664p.

1. sgm52a. plane strain tri3. check DOFs
2. sgm52b. plane strain tri3. run
3. sgm57. plane strain tri15. qn given
4. sgm511. plane strain qua4. disp given
5. sgm515. plane strain qua8. qn given
6. sgm527. plane strain qua9. qn given
7. sgm517. axisymmetric qua4
8. sgm524. 3d hex8
9. sgm530. 3d tet4
10. sgm422. small 3D frame

## Smooth contact element

1. contact01b

## General Solid Element
1. sigini01. zero displacements. initial stresses
2. sigini02. initial stresses. run simulation
3. square01. ini stress free square
4. selfweight01. self-weight
5. selfweight02. self-weight

## De Souza Neto, Peric and Owen's Book

*Reference*

De Souza Neto EA, Peric D, Owen DRJ (2008) Computational Methods For Plasticity, Wiley, 791p

1. spo751a. cylinder expansion. check DOFs
2. spo751b. cylinder expansion. run
3. spo751re. cylin exp. Richardson extrapolation
