#!/bin/bash
#
# Given a bunch of fepout files named like foobar000.fepout, calculate the total delta_G
# by invoking "assemble_and_calc_deltag.sh foobar". Assumes your mmdevel directory is in
# $HOME/mmdevel
thename=$1

MMDEVEL=$HOME/mmdevel

if [ ! -d $MMDEVEL ]
then
	echo "I thought your mmdevel directory was in $MMDEVEL. Guess not?"
	exit 1
fi

if [ -z $thename ]
then
	echo "Please specify the base name of the fepout files you are interested in."
	exit 1
fi

$MMDEVEL/fepout_cat.py --generate-spec ${thename}???.fepout > fepcalc_${thename}.yaml
$MMDEVEL/fepout_cat.py --spec fepcalc_${thename}.yaml *.fepout > all.${thename}.fepout
$MMDEVEL/calc_fepout_deltag.sh all.${thename}.fepout
