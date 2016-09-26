#!/bin/bash

FILES="upp_3mcol_desic upp_3mcol_desic16e"

for f in $FILES; do
    mpirun -np 4 gofem $f
    go run doplot.go $f
done
