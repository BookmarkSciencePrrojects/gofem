#!/bin/bash

FILES="*.go *.sim"

while true; do
    inotifywait -q -e modify $FILES
    echo
    echo
    gofem linear-qua9co && go run doplot.go
done
