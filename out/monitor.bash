#!/bin/bash

FILES="out.go results.go beamdiagrams.go t_beamdiagrams_test.go"

while true; do
    inotifywait -q -e modify $FILES
    echo
    echo
    go test -test.run="beamdiag01"
done
