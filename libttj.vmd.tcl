# Load all frames of a particular trajectory prefix
proc load_all_frames { {prefix "prod"} } {
    foreach i [lsort -dictionary [glob "$prefix*.dcd"]] {
        mol addfile $i step 10 waitfor all
    }
}


# Load structure and trajectory in a given directory, assuming you used CHARMM-GUI and
# this is a membrane protein
proc load_membrane_protein_system {dirname {prefix "prod"}} {
    set oldcwd [pwd]
    cd "$dirname"
    mol new step5_assembly.xplor_ext.psf
    mol rename top "$dirname"
    mol delrep 0 top
    mol addrep top
    mol modstyle 0 top "NewCartoon"
    mol modselect 0 top "protein"
    mol addrep top
    mol modstyle 1 top "Lines"
    mol modselect 1 top "not (protein or lipid or resname CHL1 or water or ion)"
    cd "namd"
    load_all_frames "$prefix"
    cd "$oldcwd"
}

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

# Rotates a group of atoms (in the form of an atomselect), useful for
# positioning a ligand in its site
proc rotateatoms {m x y z} {
  set loc [measure center $m]
  $m moveby [vecinvert $loc]
  $m move [transaxis x $x]
  $m move [transaxis y $y]
  $m move [transaxis z $z]
  $m moveby $loc
}

# Draw a box.
# Useful for visualizing Autodock Vina bounding boxes
proc draw_box {center_x center_y center_z size_x size_y size_z {color "yellow"}} {
    set min_x [expr $center_x - ($size_x/2)]
    set max_x [expr $center_x + ($size_x/2)]
    set min_y [expr $center_y - ($size_y/2)]
    set max_y [expr $center_y + ($size_y/2)]
    set min_z [expr $center_z - ($size_z/2)]
    set max_z [expr $center_z + ($size_z/2)]

    # Lifted from: https://www.ks.uiuc.edu/Research/vmd/vmd-1.7.1/ug/node173.html
    draw materials off
    draw color $color
    draw line "$min_x $min_y $min_z" "$max_x $min_y $min_z"
    draw line "$min_x $min_y $min_z" "$min_x $max_y $min_z"
    draw line "$min_x $min_y $min_z" "$min_x $min_y $max_z"
    draw line "$max_x $min_y $min_z" "$max_x $max_y $min_z"
    draw line "$max_x $min_y $min_z" "$max_x $min_y $max_z"
    draw line "$min_x $max_y $min_z" "$max_x $max_y $min_z"
    draw line "$min_x $max_y $min_z" "$min_x $max_y $max_z"
    draw line "$min_x $min_y $max_z" "$max_x $min_y $max_z"
    draw line "$min_x $min_y $max_z" "$min_x $max_y $max_z"
    draw line "$max_x $max_y $max_z" "$max_x $max_y $min_z"
    draw line "$max_x $max_y $max_z" "$min_x $max_y $max_z"
    draw line "$max_x $max_y $max_z" "$max_x $min_y $max_z"
}