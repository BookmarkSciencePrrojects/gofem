#!/bin/bash

FILES="msh.go sim.go t_nurbs_test.go t_read_test.go"

while true; do
    inotifywait -q -e modify $FILES
    echo
    echo
    go test -test.run="msh02"
done
