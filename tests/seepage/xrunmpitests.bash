#!/bin/bash

go build -o /tmp/gofem/liquid_main liquid_main.go && mpirun -np 3 /tmp/gofem/liquid_main
