#!/bin/bash

set -e

ALL="
dynamics_sgbook
patch_test
rjoint_ex01_curved
rjoint_ex06_pullout
seep_ex01_freesurf
seep_ex02_freesurf
seep_simple_flux
spo751_pressurised_cylinder
spo754_strip_footing_collapse
up_3mcolumn_desiccation
up_indentation2d_unsat
upp_3mcolumn_desiccation
"

# sed -i -e 's/utl.Panic/chk.Panic/g' \
#       -e 's/utl.Err\>/chk.Err/g' \
#       -e 's/utl.Sramp/fun.Sramp/g' $f

runcommand() {
    exe=$1
    echo
    echo
    echo ">>>>>>>>>>>>>>>> $exe <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    for f in *.go; do
        if [[ -f $f ]]; then
            echo $f
            sed -i 's/fun.Prm/fun.P/g' $f
            goimports -w $f
        fi
    done
}

for exe in $ALL; do
    HERE=`pwd`
    cd "examples/"$exe
    runcommand $exe
    cd $HERE
done
