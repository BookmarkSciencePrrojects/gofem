#!/bin/bash

set -e

echo "usage:"
echo "    $0 JOB"
echo "where JOB is:"
echo "    0 -- count lines [default]"
echo "    1 -- execute goimports"
echo "    2 -- generate depedency graphs"
echo "    3 -- fix links in README files"
echo "    4 -- run sed [must be customised]"

JOB=0
if [[ $# != 0 ]]; then
    JOB=$1
    if [[ $JOB -lt 0 || $JOB -gt 4 ]]; then
        echo
        echo "Job number $1 is invalid"
        echo
        exit 1
    fi
fi

echo "current JOB = $JOB"

if [[ $JOB == 0 ]]; then
    totnfiles=0
    totnlines=0
    for f in `find . -iname "*.go"`; do
        totnfiles=$(($totnfiles+1))
        totnlines=$(($totnlines+`wc -l $f | awk '{print $1}'`))
    done
    echo
    echo "Total number of files = $totnfiles"
    echo "Total number of lines = $totnlines"
    exit 0
fi

ALL=" \
ana \
shp \
mdl/generic \
mdl/solid \
mdl/fluid \
mdl/conduct \
mdl/retention \
mdl/diffusion \
mdl/thermomech \
mdl/porous \
inp \
ele \
ele/solid \
ele/seepage \
ele/diffusion \
ele/thermomech \
ele/porous \
fem \
tests \
tests/solid \
tests/seepage \
tests/diffusion \
tests/thermomech \
tests/porous \
out \
"

EXA=" \
examples/dynamics_sgbook \
examples/patch_test \
examples/rjoint_ex01_curved \
examples/rjoint_ex06_pullout \
examples/seep_ex01_freesurf \
examples/seep_ex02_freesurf \
examples/seep_simple_flux \
examples/spo751_pressurised_cylinder \
examples/spo754_strip_footing_collapse \
examples/up_3mcolumn_desiccation \
examples/up_indentation2d_unsat \
examples/upp_3mcolumn_desiccation \
"

rungoimports() {
    pkg=$1
    for f in *.go; do
        echo $f
        goimports -w $f
    done
}

depgraph(){
    pkg=$1
    fna="/tmp/gofem/depgraph-${pkg/\//_}-A.png"
    fnb="/tmp/gofem/depgraph-${pkg/\//_}-B.svg"
    godepgraph -s github.com/cpmech/gofem/$pkg | dot -Tpng -o $fna
    graphpkg -stdout -match 'gofem' github.com/cpmech/gofem/$pkg > $fnb
    echo "file <$fna> generated"
    echo "file <$fnb> generated"
}

fixreadme() {
    pkg=$1
    old="http://rawgit.com/cpmech/gofem/master/doc/xx${pkg/\//-}.html"
    new="https://godoc.org/github.com/cpmech/gofem/${pkg}"
    sed -i 's,'"$old"','"$new"',' README.md
}

custom() {
    pkg=$1
    old="fun.Params"
    new="dbf.Params"
    for f in *.go; do
        echo $f
        sed -i 's,'"$old"','"$new"',' $f
        goimports -w $f
    done
}

if [[ $JOB == 1 ]]; then
    ALL="$ALL $EXA"
fi

if [[ $JOB == 2 ]]; then
    mkdir -p /tmp/gofem
fi

idx=1
for pkg in $ALL; do
    HERE=`pwd`
    cd $pkg
    echo
    echo ">>> $idx $pkg <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    if [[ $JOB == 1 ]]; then
        rungoimports $pkg
    fi
    if [[ $JOB == 2 ]]; then
        depgraph $pkg
    fi
    if [[ $JOB == 3 ]]; then
        fixreadme $pkg
    fi
    if [[ $JOB == 4 ]]; then
        custom $pkg
    fi
    cd $HERE
    (( idx++ ))
done
