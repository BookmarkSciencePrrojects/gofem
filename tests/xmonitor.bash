#!/bin/bash

FILES="../fem/*.go ../ele/*.go *.go"

dorun(){
    echo
    echo
    go test
}

while true; do
    inotifywait -q -e modify $FILES
    dorun
done
