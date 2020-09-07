# Load all frames of a particular trajectory prefix
proc load_all_frames { {prefix "prod"} {step 10} } {
    set files [lsort -dictionary [glob -nocomplain "$prefix*.dcd"]] 
    if {[llength $files] == 0} {
        puts "No DCD files with prefix $prefix were found."
    }
    foreach i $files {
        # Use a regexp to avoid catching things like "prodB1" when the user specified only "prod"
        if {[regexp "$prefix\[0-9\]+" $i]} {
            mol addfile $i step $step waitfor all
        } else {
            puts "load_all_frames: Not loading $i since it isn't $prefix and a number."
        }
    }
}


# Load structure and trajectory in a given directory, assuming you used CHARMM-GUI and
# this is a membrane protein
proc load_membrane_protein_system {dirname {prefix "prod"} {step 10}} {
    set oldcwd [pwd]
    cd "$dirname"
    # These days CHARMM-GUI names the PSF file differently, so check for new and old
    set psf [glob "namd/step5_charmm2namd.psf" "step5_assembly.xplor_ext.psf"]
    mol new "$psf"
    mol rename top "${dirname}_${prefix}"
    mol delrep 0 top
    mol addrep top
    mol modstyle 0 top "NewCartoon"
    mol modselect 0 top "protein"
    mol addrep top
    mol modstyle 1 top "Lines"
    mol modselect 1 top "not (protein or lipid or resname CHL1 or water or ion)"
    cd "namd"
    load_all_frames "$prefix" $step
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

# Calculate per-residue RMSD over a trajectory.
# Presumes that proteins are the same length.
# The ref_molid frame is not advanced - only that of molid
proc per_residue_rmsd {ref_molid ref_segid ref_residues molid {outfile "-"}} {
    set all_ref_residues [lsort -unique -dictionary [[atomselect $ref_molid "protein"] get resid]]
    set all_target_residues [lsort -unique -dictionary [[atomselect $molid "protein"] get resid]]
    set num_target_residues [llength $all_target_residues]
    set num_ref_residues [llength $all_ref_residues]

    if {$num_ref_residues != $num_target_residues} {
        puts "Sequence lengths of molids $ref_molid ($num_ref_residues) and $molid ($num_target_residues) don't match."
        puts "Which means they might not be the same protein."
        puts "I'm not smart enough to decide which residues should be compared."
        return 1
    }

    # Provide shorthand to specify alignment on all residues
    if {$ref_residues == "all"} {
        set ref_residues $all_ref_residues
    }

    # Cache atomselect objects for all residues in both proteins
    set column_labels {}
    set all_segids [lsort -unique [[atomselect $ref_molid protein] get segid]]
    set ref_res_sel {}
    set ws_res_sel {}
    foreach segid $all_segids {
        foreach resid $all_ref_residues {
            set ref_res [atomselect $ref_molid "protein and segid $segid and resid $resid"]
            set ws_res [atomselect $molid "protein and segid $segid and resid $resid"]
            if {[llength [$ref_res list]] == 0} {
                puts "Info: ref $segid:$resid does not exist"
                continue
            }
            if {[llength [$ws_res list]] == 0} {
                puts "Info: target $segid:$resid does not exist"
                continue
            }
            lappend ref_res_sel $ref_res
            lappend ws_res_sel $ws_res
            # Construct labels for columns in the output CSV file
            set ref_res_name [string totitle [lindex [$ref_res get resname] 0]]
            lappend column_labels "$ref_res_name${resid}_${segid}"
        }
    }

    # We'll dump results into a CSV file
    if {$outfile == "-"} {
        set outfile "rmsd_${ref_molid}_vs_${molid}.csv"
    }
    set fd [open $outfile "w"]

    # Print CSV column headers
    lappend column_labels "AllFitResidues"
    puts $fd [join $column_labels ","]

    # Select the atoms we're going to use for the fit
    # We do this outside the loop because each time you do this, a new object is created
    set ref_fit [atomselect $ref_molid "protein and segid $ref_segid and name CA and resid $ref_residues"]
    # set ref_all [atomselect $ref_molid all]
    set target_fit [atomselect $molid "protein and segid $ref_segid and name CA and resid $ref_residues"]
    set target_all [atomselect $molid all]

    puts "Number of reference residues used in fit: [llength [$ref_fit get name]]"
    puts "Number of target residues used in fit: [llength [$target_fit get name]]"
    
    # Iterate over frames of target. On each frame do a fit on $target and $ref
    # against whatever frame the reference is on
    set num_res_sel [llength $ref_res_sel]
    set num_frames [molinfo $molid get numframes]
    for {set frame 0} {$frame < $num_frames} {incr frame} {
        # Set the current frame
        $target_all frame $frame
        $target_fit frame $frame
        # Do the alignment
        set trans_mat [measure fit $target_fit $ref_fit]
        $target_all move $trans_mat
        set all_fit_rmsd [measure rmsd $target_fit $ref_fit]
        set per_res_rmsd {}
        for {set i 0} {$i < $num_res_sel} {incr i} {
            [lindex $ws_res_sel $i] frame $frame
            lappend per_res_rmsd [format "%.2f" [measure rmsd [lindex $ref_res_sel $i] [lindex $ws_res_sel $i]]]
        }
        # Print the individual residue RMSDs
        puts -nonewline $fd [join $per_res_rmsd ","]
        # Print the overall RMSD of all residues in question
        puts -nonewline $fd ","
        puts $fd "$all_fit_rmsd"
    }

    # Close the output file
    close $fd
    puts "Wrote output per-frame RMSD in CSV format to: $outfile"

    # Delete all those atomselect objects because apparently there is no garbage collection
    # $ref delete
    $ref_fit delete
    $target_all delete
    $target_fit delete
    foreach sel $ref_res_sel { $sel delete }
    foreach sel $ws_res_sel { $sel delete }
}


# Tcl does not include a range function. In this day and age.
proc range {from to} {
    set out {}
    for {set i $from} {$i < $to} {incr i} {
        lappend out $i
    }
    return $out
}


# Automate part of the process of aligning one trajectory to another, when
# the proteins are similar but not structurally exactly the same.
# In VMD, the easiest way I can see to do this normally is using the
# STAMP structural alignment feature of the MultiSeq plugin to rotate
# the first frame of the to-be-aligned trajectory, then align that
# trajectory against that structure.
# This function automates what STAMP is used for here, but not using
# STAMP. Instead we use UCSF Chimera MatchMaker.
proc generate_struct_align_ref {to_align_molid ref_molid} {
    animate goto start
    [atomselect $ref_molid protein] writepdb ref_prot.tmp.pdb
    [atomselect $to_align_molid protein] writepdb to_align_prot.tmp.pdb
    set f [open "align.tmp.cmd" "w"]
    puts $f "mmaker #0 #1"
    puts $f "write relative #0 #1 align_molid${to_align_molid}_to_me.pdb"
    close $f
    exec chimera --nogui --silent ref_prot.tmp.pdb to_align_prot.tmp.pdb align.tmp.cmd
    mol new "align_molid${to_align_molid}_to_me.pdb"
    file delete ref_prot.tmp.pdb to_align_prot.tmp.pdb align_tmp.cmd
    file delete "align_molid${to_align_molid}_to_me.pdb"
    puts "Generated and loaded a proxy reference structure for molid ${to_align_molid}."
}


# Substitute a canonical name for the various states of histidine,
# so that these don't cause align_trajectory to think the proteins are different
proc sanitize_resname_list {resnames} {
    set out_list ""
    foreach res $resnames {
        if {$res == "HSD" || $res == "HSP" || $res == "HSE"} {
            lappend out_list "HIS"
        } else {
            lappend out_list $res
        }
    }
    return $out_list
}


# Aligns a trajectory so as to remove global rotations and translations
# Will use the first frame of ref_molid as reference. If not specified, use the
# first frame of the trajectory
proc align_trajectory {to_align_molid {ref_molid "same"} {align_sel "backbone"}} {
    if {$ref_molid eq "same"} {
        set ref_molid $to_align_molid
    }

    puts "Aligning molid $to_align_molid against the first frame of reference molid $ref_molid."

    set ref_sel [atomselect $ref_molid $align_sel]
    set to_align_sel [atomselect $to_align_molid $align_sel]
    set proteins_are_same 1

    # Are these actually the same protein?
    # If not, generate a temporary proxy reference structure that we can use "measure fit" with
    if {[sanitize_resname_list [$ref_sel get resname]] != [sanitize_resname_list [$to_align_sel get resname]]} {
        puts "But these are not the same protein because they have different resnames."
        puts "Therefore, generating a proxy reference structure by structural alignment using Chimera MatchMaker."
        generate_struct_align_ref $to_align_molid $ref_molid
        # No need for the old reference atomselect...we need to use the new one now
        $ref_sel delete
        set ref_sel [atomselect top $align_sel]
        set proteins_are_same 0
    }

    # Specify which atoms to move to best fit position: all of them
    set move_sel [atomselect $to_align_molid "all"]
    # Use only the first frame of the reference trajectory
    $ref_sel frame 0

    # Iterate over frames of to_align_molid
    set num_frames [molinfo $to_align_molid get numframes]
    for {set frame 0} {$frame < $num_frames} {incr frame} {
        # Select the frame of the trajectory we will be fitting
        $to_align_sel frame $frame
        $move_sel frame $frame
        # Calculate the transformation matrix for the fit
        set mat [measure fit $to_align_sel $ref_sel]
        # Actually move the atoms
        $move_sel move $mat
    }
    # Clean up our atomselects, because we are a good citizen
    $ref_sel delete
    $to_align_sel delete
    $move_sel delete
    # If we generated a temporary proxy reference structure, delete it
    if {$proteins_are_same eq 0} {
        mol delete top
    }
    puts "Done with trajectory alignment. Enjoy!"
}

# Computes the RMSD of something over time from frame 0
# Does NOT do a fit because maybe you wanted a global alignment done a particular way.
# Or some other weird thing. So the caller needs to do that first.
proc rmsd_over_time_nofit { {mol top} {seltext "noh"} } {
    set ref [atomselect $mol "$seltext" frame 0]
    set thing [atomselect $mol "$seltext"]
    set num_frames [molinfo $mol get numframes]
    set rmsd {}
    for {set i 0} {$i < $num_frames} {incr i} {
        $thing frame $i
        lappend rmsd [measure rmsd $ref $thing]
    }
    # Be a good citizen and delete our atomselect objects
    $ref delete
    $thing delete
    return $rmsd
}

# Write a list to a single column CSV file
# Although with only one column there are no commas
proc write_to_csv {data headername filename} {
    # We'll dump results into a CSV file
    set fd [open $filename "w"]

    # Write CSV column header
    puts $fd $headername

    # Write the data
    puts $fd [join $data "\n"]
    close $fd
}
