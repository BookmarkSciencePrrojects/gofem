#!/bin/bash

FILES="*.go"

while true; do
    inotifywait -q -e modify $FILES
    echo
    echo
    echo
    echo
    go test -test.run="hyperelast03"
done
