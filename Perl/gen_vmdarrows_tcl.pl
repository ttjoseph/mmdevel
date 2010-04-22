#!/usr/bin/env perl
#
# Generates a TCL script to be run by VMD that draws displacement arrows
# from essential dynamics (or something else.)
#
# Format of each input file is three real numbers per line (that is, XYZ).
#
# Usage: <structure> <arrow-vectors>

open STRUCT, $ARGV[0] or die "Couldn't open structure file $ARGV[0]\n";
open VECS, $ARGV[1] or die "Couldn't open vectors file $ARGV[1]\n";

$scale = 70;

print "proc vmd_draw_arrow {mol start end} {\n";
print "   # Shamelessly ripped off from a FORTRAN program by E. Giudice\n";
print "   # An arrow is made of a cylinder and a cone\n";
print "   set middle [vecadd \$start [vecscale 0.8 [vecsub \$end \$start]]]\n";
print "   graphics \$mol cylinder \$start \$middle radius 0.15\n";
print "   graphics \$mol cone \$middle \$end radius 0.25\n";
print "}\n\n";

while(<STRUCT>)
{
	s/^\s+//g; # Get rid of stray whitespace
	s/\s+/ /g;
	($index, $x1, $y1, $z1) = split /\s+/;
	#print STDERR "$x1 $y1 $z1\n";
	$_ = <VECS>;
	s/^\s+//g;
	s/\s+/ /g;
	($dx, $dy, $dz) = split /\s+/;
	#print STDERR "...$dx///$dy---$dz\n";
	$x2 = $x1 + $scale * $dx;
	$y2 = $y1 + $scale * $dy;
	$z2 = $z1 + $scale * $dz;

	# Function naming in TCL and VMD seems pretty bizarre.
	# I thought it was called vmd_draw_arrow?
	print "draw arrow { $x1 $y1 $z1 } { $x2 $y2 $z2 }\n";
}
