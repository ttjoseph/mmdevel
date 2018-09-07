#!/usr/bin/bash
#
# Wraps a trajectory back into one PBC cell, using qwrap.
#
# Writes the result to wrapped.dcd because I'm too lazy to do proper
# command line arguments. This script is just to save little bit of
# typing.

psf=$1
shift
dcds=$*

echo === Topology file ===
echo $psf
echo === Trajectory files ===
echo $dcds

vmd -dispdev text <<HOORAY
package require qwrap
mol new $psf
foreach i {$dcds} { mol addfile \$i waitfor all }
qwrap
[atomselect top all] writedcd wrapped.dcd
quit
HOORAY

echo === Wrapped trajectory written to wrapped.dcd ===