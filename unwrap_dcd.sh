#!/bin/bash

CPPTRAJ=$HOME/build/amber16/bin/cpptraj
psf=$1
shift

$CPPTRAJ <<FOO
parm $psf
$(echo $* | tr ' ' '\n' | awk '{print "trajin " $1}')
unwrap
trajout $(basename $psf .psf).unwrapped.dcd
go
FOO
