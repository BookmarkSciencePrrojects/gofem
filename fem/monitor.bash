#!/bin/bash

FILES="*.go data/*.*"

while true; do
    inotifywait -q -e modify $FILES
    echo
    echo
    go test -test.run="upp01a"
done
