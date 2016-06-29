#!/bin/bash

FILES="*.go"

dorun(){
    echo
    echo
    echo
    echo
    go test -test.run="upp01"
}

while true; do
    inotifywait -q -e modify $FILES
    dorun
done
