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

compile ana
compile shp
compile mdl/sld
compile mdl/cnd
compile mdl/lrm
compile mdl/fld
compile mdl/por
compile inp
compile fem
compile out

echo
echo
echo "[1;32m>>> compiling binaries <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<[0m"
make
