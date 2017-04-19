#!/bin/bash

set -e

GEN="ana ele fem inp out shp tests tools"
ELE="ele/diffusion ele/porous ele/seepage ele/solid ele/thermomech"
MDL="mdl/conduct mdl/diffusion mdl/fluid mdl/generic mdl/porous mdl/retention mdl/solid mdl/thermomech"
TST="tests/diffusion tests/porous tests/seepage tests/solid tests/thermomech"
EXA="examples/dynamics_sgbook examples/patch_test examples/rjoint_ex01_curved examples/rjoint_ex06_pullout examples/seep_simple_flux examples/spo751_pressurised_cylinder examples/spo754_strip_footing_collapse examples/up_3mcolumn_desiccation examples/up_indentation2d_unsat examples/upp_3mcolumn_desiccation"
ALL="$GEN $ELE $MDL $TST $EXA"

runcommand() {
    pkg=$1
    echo
    echo
    echo ">>>>>>>>>>>>>>>> $pkg <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    for f in *.go; do
        echo $f
        sed -i -e '/SetForPng/d' $f
        sed -i -e '/SetForEps/d' $f
        goimports -w $f
    done
}

for pkg in $ALL; do
    HERE=`pwd`
    cd $pkg
    runcommand $pkg
    cd $HERE
done
