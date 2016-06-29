#!/bin/bash

FILES="../../fem/*.go ../../ele/*.go ../*.go *.go"

dorun(){
    echo
    echo
    echo
    echo
    go test
    #go test -test.run="upp01a"
}

while true; do
    inotifywait -q -e modify $FILES
    dorun
done
