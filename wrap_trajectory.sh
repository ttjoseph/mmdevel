#!/usr/bin/bash
#
# Wraps a trajectory back into one PBC cell, using qwrap.
# If you specify more than one trajectory, each one is written to its
# own wrapped trajectory.

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
foreach i {$dcds} {
	animate delete all
	mol addfile \$i waitfor all
	set basename [file rootname \$i]
	set wrappedname "\${basename}.wrapped.dcd"
	qwrap
	animate write dcd \$wrappedname
}
quit
HOORAY

