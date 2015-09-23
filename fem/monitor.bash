#!/bin/bash

FILES="e_beam.go e_u.go e_p.go e_up.go t_nurbs_test.go t_beams_test.go"

while true; do
    inotifywait -q -e modify $FILES
    echo
    echo
    go test -test.run="beam04"
done
