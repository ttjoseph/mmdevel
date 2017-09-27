# Get VMD to align two local ligand environments to each other, despite
# its best efforts not to be smart enough to do that.
#
# The motivating use case is aligning two of the same protein subunit together that have
# different ligands bound (or one with and one without a ligand).
#
# Example:
#     set_beta_for_alignment 0 "segname PROA and resid 100 105 136" 1 "segname PROB and resid 100 105 136"
#     (Then use RMSD Trajectory Tool to align on "beta 1")

proc set_beta_for_alignment {mol1 spec1 mol2 spec2} {
 	[atomselect $mol1 all] set beta 0
 	[atomselect $mol2 all] set beta 0
 	set mol1_sel [atomselect $mol1 $spec1]
 	set mol2_sel [atomselect $mol2 $spec2]
 	if {[llength $mol1_sel] == [llength $mol2_sel] } {
	 	$mol1_sel set beta 1
	 	$mol2_sel set beta 1
	 	puts "OK, now you can align those things on 'beta 1' using RMSD Trajectory Tool"
 	} else {
 		puts "Those selections have different numbers of atoms and RMSD Trajectory Tool is too dumb for that."
 	}

} 