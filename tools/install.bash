#!/bin/bash

GP="${GOPATH%:*}"

echo "...installing tools to $GP"

go build -o /tmp/gofem/GenVtu GenVtu.go && mv /tmp/gofem/GenVtu $GP/bin/
echo "......GenVtu installed"

go build -o /tmp/gofem/MatTable MatTable.go && mv /tmp/gofem/MatTable $GP/bin/
echo "......MatTable installed"

go build -o /tmp/gofem/PlotLrm PlotLrm.go && mv /tmp/gofem/PlotLrm $GP/bin/
echo "......PlotLrm intalled"

go build -o /tmp/gofem/LocCmDriver LocCmDriver.go && mv /tmp/gofem/LocCmDriver $GP/bin/
echo "......LocCmDriver installed"

go build -o /tmp/gofem/ResidPlot ResidPlot.go && mv /tmp/gofem/ResidPlot $GP/bin/
echo "......ResidPlot installed"

go build -o /tmp/gofem/Msh2vtu Msh2vtu.go && mv /tmp/gofem/Msh2vtu $GP/bin/
echo "......Msh2vtu installed"
