#!/bin/bash

FILES="*.tex"

while true; do
    inotifywait -q -e modify $FILES
    echo
    echo
    make
done
