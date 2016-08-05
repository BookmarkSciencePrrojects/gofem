#!/bin/bash

set -e

HERE=`pwd`

testandinstall() {
    echo
    echo
    echo "[1;32m>>> test-and-install $1 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<[0m"
    cd $1
    touch *.go
    go test
    go install
    cd $HERE
}

testonly() {
    echo
    echo
    echo "[1;32m>>> test-only $1 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<[0m"
    cd $1
    touch *.go
    go test
    cd $HERE
}

installonly() {
    echo
    echo
    echo "[1;32m>>> install-only $1 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<[0m"
    cd $1
    touch *.go
    go install
    cd $HERE
}


testandinstall ana
testandinstall shp
testandinstall mdl/generic
testandinstall mdl/solid
testandinstall mdl/fluid
testandinstall mdl/conduct
testandinstall mdl/retention
testandinstall mdl/diffusion
testandinstall mdl/thermomech
testandinstall mdl/porous
testandinstall inp
testandinstall ele
testandinstall ele/solid
testandinstall ele/seepage
testandinstall ele/diffusion
#testandinstall ele/thermomech
testandinstall ele/porous
testandinstall fem
testandinstall tests
testonly       tests/solid
testonly       tests/seepage
testonly       tests/diffusion
testonly       tests/thermomech
testonly       tests/porous
testandinstall out

echo
echo
echo "[1;32m>>> compiling binaries <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<[0m"
make
