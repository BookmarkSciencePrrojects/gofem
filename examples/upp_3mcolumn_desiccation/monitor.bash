#!/bin/bash

FILES="*.go *.sim"

while true; do
    inotifywait -q -e modify $FILES
    echo
    echo
    gofem upp_3mcol_desic && go run doplot.go
done
