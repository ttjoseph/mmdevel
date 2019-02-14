#!/bin/sh

FEP_BASENAME=$1
shift

python ~/mmdevel/fepout_cat.py ${FEP_BASENAME}???.fepout > $FEP_BASENAME.all.fepout
grep FepEnergy $FEP_BASENAME.all.fepout | awk '{print $10}' > $FEP_BASENAME.fepenergy.gnuplot
# 8 12
cat ${FEP_BASENAME}*.fepout | grep '#Free' | awk '{print $8" "$12}' > $FEP_BASENAME.delta.gnuplot

gnuplot <<FOO
set terminal svg
set output "$FEP_BASENAME.fepenergy.svg"
plot "$FEP_BASENAME.fepenergy.gnuplot" w l

set output "$FEP_BASENAME.delta.svg"
plot "$FEP_BASENAME.delta.gnuplot" smooth cumulative
FOO
