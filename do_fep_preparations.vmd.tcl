# Helper stuff for generating the various files necessary for a FEP with spherical restraints

# Assumes the top molecule is the system and we want to disappear $ligspec
# and that the ligand is not a protein
proc do_fep_preparations {ligspec} {
	# Make disappear.fep
	set all [atomselect top all]
	set lig [atomselect top "$ligspec"]

	$all set beta 0
	$all set occupancy 0

	# Tell NAMD which atoms to decouple
	$lig set beta -1
	$all writepdb disappear.fep
	puts "# Wrote FEP atom selection disappear.fep"

	# Tell colvars which atoms to use for reorientation
	$all set beta 0
	set restref [atomselect top "protein and backbone within 7 of ($ligspec)"]
	set restref_len [llength [$restref list]]
	puts "# rest_ref uses $restref_len atoms. Make sure that is OK!"
	puts ""

	$restref set occupancy 1
	$all writepdb rest_ref.pdb
	puts "# Wrote restraint reference rest_ref.pdb"
	$all set occupancy 0
	[atomselect top "protein and name CA"] set occupancy 1
	$all writepdb dont_spin_ref.pdb
	puts "# Wrote dont_spin_ref.pdb"

	# Show ligand atomNumbers for the benefit of colvars
	set nums {}
	foreach x [$lig list] {
		lappend nums [expr {$x + 1}]
	}
	puts "atomNumbers $nums"

	# Show COM (COG?) of the ligand, again for colvars
	puts "# Center of ligand $ligspec:"
	set lig_center [join [measure center $lig] ", "]
	puts "dummyAtom ($lig_center)"
}