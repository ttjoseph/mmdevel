#!/bin/bash
#
# Convenience script to plot the output of angle_over_traj.py,
# assuming there are two columns

gnuplot -persist <<FOO
set key autotitle columnhead
set datafile separator ','
set xlabel 'Frame (1000 frames per 20 ns)'
set ylabel 'Angle (degrees)'
plot '$1' u 0:1 w l, '' u 0:2 w l
FOO
