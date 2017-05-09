#!/bin/bash

vmd -dispdev text <<FOO
package require parsefep
parsefep -gc 0 -forward $1
quit
FOO

grep FepEnergy $1 | awk '{print $10}' > $1.fepenergy.gnuplot

echo "Generated Gnuplot plot in $1.fepenergy.gnuplot"
