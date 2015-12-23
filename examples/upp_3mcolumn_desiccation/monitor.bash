#!/bin/bash

FILES="*.go *.sim *.mat *.msh"

while true; do
    inotifywait -q -e modify $FILES
    echo
    echo
    gofem upp_3mcol_desic16e && go run doplot.go
done
