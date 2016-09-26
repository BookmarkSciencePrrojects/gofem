#!/bin/bash

rm -rf /tmp/gofem

nerrors=0

for i in {1..100}; do
    echo $i
    mpirun -np 4 gofem upp_3mcol_desic false
    res=`file /tmp/gofem/upp_3mcol_desic/upp_3mcol_desic_p0_nod_0000000000.gob | grep 'cannot open'`
    if [[ ! -z "$res" ]]; then
        echo "error"
        nerrors=$((nerrors+1))
    fi
done

echo
echo
echo "number of errors = ${nerrors}"
