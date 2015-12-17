#!/bin/bash

FILES="*.go *.sim *.msh *.mat"

while true; do
    inotifywait -q -e modify $FILES
    echo
    echo
    gofem liqgas-3mcol-4el.sim
    go run doplot.go
done
