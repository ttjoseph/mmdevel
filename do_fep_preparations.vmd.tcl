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
	# Try to use a little under 80 atoms for the fit
	$all set beta 0
	set restref_len 1000
	set restref_radius 8
	# The for loop construction in Tcl appears to be brain-damaged, so we use a while loop
	while {$restref_len > 80} {
		set restref [atomselect top "protein and backbone within $restref_radius of ($ligspec)"]
		set restref_len [llength [$restref list]]
		set restref_radius [expr {$restref_radius - 0.1}]
	}

	puts "# rest_ref uses $restref_len atoms. Make sure that is OK!"
	puts ""

	$restref set occupancy 1
	$all writepdb rest_ref.pdb
	puts "# Wrote restraint reference rest_ref.pdb"
	$all set occupancy 0

	# Reference structure to keep the protein from spinning
	[atomselect top "protein and name CA"] set occupancy 1
	$all writepdb dont_spin_ref.pdb
	puts "# Wrote dont_spin_ref.pdb"

	# Show ligand atomNumbers for the benefit of colvars
	set nums {}
	foreach x [$lig list] {
		lappend nums [expr {$x + 1}]
	}
	# Show COM (COG?) of the ligand, again for colvars
	set lig_center [join [measure center $lig] ", "]

	set fd [open "restraints.ligand_fb.col" w]
	puts $fd "# Colvars restraint setup for FEP for ligand \"$ligspec\"

# Distance of ligand from its binding site
colvar {
    name ligand_fb

    width 1.0
    upperBoundary 5.0
    upperWallConstant 100 

    distance {
        group1 {
            # Ligand atom numbers
            atomNumbers $nums
            rotateReference
            centerReference
            # fit to $restref_len protein atoms within $restref_radius A
            fittingGroup {
                atomsCol O
                atomsFile rest_ref.pdb
            } 
            refPositionsFile rest_ref.pdb       # ref protein coords for fitting. 
        }
        group2 {
            dummyAtom ($lig_center)
        }
    }

}

# With a lambda schedule, we could fade away this restraint 
harmonic {
    colvars ligand_fb
    forceConstant 0
    centers 0
}

# Aped from https://github.com/colvars/colvars/blob/master/examples/03_orientation.colvars.in
colvar {
    name dont_spin_protein
    orientation {
        atoms {
            atomsCol O
            atomsFile dont_spin_ref.pdb
        }
        refPositionsFile dont_spin_ref.pdb 
    }
}

harmonic {
    colvars dont_spin_protein
    # Restraint center is the unit quaternion representing no rotation
    centers (1.0, 0.0, 0.0, 0.0)
    # Chosen at random more or less
    forceConstant 5.0
}
"
	puts "# Wrote restraints colvars file restraints.ligand_fb.col"
	close $fd
}