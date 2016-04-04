#!/bin/bash

FILES="*.go data/*.*"

while true; do
    inotifywait -q -e modify $FILES
    echo
    echo
    go test -test.run="diffu01b"
done
