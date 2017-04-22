#!/bin/bash

FILES="*.go data/*.sim"

while true; do
    inotifywait -q -e modify $FILES
    echo
    echo
    go test -test.run="msh04"
done
