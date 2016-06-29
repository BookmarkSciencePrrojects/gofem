#!/bin/bash

HERE=`pwd`

compile() {
    echo
    echo
    echo "[1;32m>>> compiling $1 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<[0m"
    cd $1
    touch *.go
    go test
    go install
    cd $HERE
}

#compile ana
#compile shp
#compile mdl/generic
#compile mdl/solid
#compile mdl/fluid
#compile mdl/conduct
#compile mdl/retention
#compile mdl/diffusion
#compile mdl/thermomech
#compile mdl/porous
#compile inp
compile ele
compile ele/solid
compile ele/seepage
compile ele/diffusion
#compile ele/thermomech
compile ele/porous
compile fem
compile tests
compile tests/solid
compile tests/seepage
compile tests/diffusion
compile tests/thermomech
compile tests/porous
#compile out

echo
echo
echo "[1;32m>>> compiling binaries <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<[0m"
make
