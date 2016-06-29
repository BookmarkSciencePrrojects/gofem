#!/bin/bash

FILES="*.go data/*.sim"

while true; do
    inotifywait -q -e modify $FILES
    echo
    echo
    go test -test.run="sim03"
done
