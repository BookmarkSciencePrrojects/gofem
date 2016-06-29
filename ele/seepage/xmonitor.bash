#!/bin/bash

FILES="*.go"

dorun(){
    echo
    echo
    echo
    echo
    go test
}

while true; do
    inotifywait -q -e modify $FILES
    dorun
done
