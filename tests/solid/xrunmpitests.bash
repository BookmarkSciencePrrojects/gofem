#!/bin/bash

go build -o /tmp/gofem/bh16_main bh16_main.go && mpirun -np 3 /tmp/gofem/bh16_main
go build -o /tmp/gofem/spo751_main spo751_main.go && mpirun -np 3 /tmp/gofem/spo751_main
