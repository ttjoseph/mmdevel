package require psfgen

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
    set psf [glob "namd/step5_input.psf" "namd/step5_charmm2namd.psf" "step5_assembly.xplor_ext.psf"]
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
proc do_fep_preparations {ligspec {filename_suffix {}}} {
	if {$filename_suffix ne {}} {
        set filename_suffix "${filename_suffix}."
    }
    
    # Make disappear.fep
	set all [atomselect top all]
	set lig [atomselect top "$ligspec"]

	$all set beta 0
	$all set occupancy 0

	# Tell NAMD which atoms to decouple
	$lig set beta -1
	$all writepdb disappear.${filename_suffix}fep
	puts "# Wrote FEP atom selection disappear.${filename_suffix}fep"

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
	$all writepdb "rest_ref.${filename_suffix}pdb"
	puts "# Wrote restraint reference rest_ref.${filename_suffix}pdb"
	$all set occupancy 0

	# Reference structure to keep the protein from spinning
	[atomselect top "protein and name CA"] set occupancy 1
	$all writepdb "dont_spin_ref.${filename_suffix}pdb"
	puts "# Wrote dont_spin_ref.${filename_suffix}pdb"

	# Show ligand atomNumbers for the benefit of colvars
	set nums {}
	foreach x [$lig list] {
		lappend nums [expr {$x + 1}]
	}
	# Show COM (COG?) of the ligand, again for colvars
	set lig_center [join [measure center $lig] ", "]

	set fd [open "restraints.${filename_suffix}colvars" w]
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
                atomsFile rest_ref.${filename_suffix}pdb
            } 
            refPositionsFile rest_ref.${filename_suffix}pdb       # ref protein coords for fitting. 
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
            atomsFile dont_spin_ref.${filename_suffix}pdb
        }
        refPositionsFile dont_spin_ref.${filename_suffix}pdb 
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
	puts "# Wrote restraints colvars file restraints.${filename_suffix}colvars"
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
proc per_residue_rmsd {ref_molid ref_segid ref_residues molid {outfile "-"} {which_atoms "name CA"}} {
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
            set ref_res [atomselect $ref_molid "protein and segid $segid and resid $resid and $which_atoms"]
            set ws_res [atomselect $molid "protein and segid $segid and resid $resid and $which_atoms"]
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
    file delete ref_prot.tmp.pdb to_align_prot.tmp.pdb align.tmp.cmd
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

# Removes duplicates without sorting the input list - adapted from wiki.tcl.tk
proc lunique {L} {
    set unique_elements {}
    foreach elt $L {
        if {[lsearch -exact $unique_elements $elt] == -1} {
            lappend unique_elements $elt
        }
    }
    return $unique_elements
}

# Reads PSF/PDB, mutates one residue, then writes a new PSF and PDB.
# Useful for relative FEP calculations among other things.
# Example:
# mutate_psf {toppar/toppar_water_ions.str top_all36_lipid_alchemy.rtf} step5_assembly.psf step5_assembly.pdb MEMB 1 POEA foo
proc mutate_psf {rtfs psf pdb mutate_segid mutate_resid mutate_to out_prefix} {
    # Split the supplied system into its constituent segments, because
    # psfgen is too stupid to deal with them all together, from what I can tell
    mol new $psf
    mol addfile $pdb waitfor all

    set segids [lunique [[atomselect top all] get segid]]
    set seg_pdbs {}

    # Load all our topologies
    resetpsf
    foreach rtf $rtfs {
        topology $rtf
    }

    foreach segid $segids {
        set seg_pdb "${out_prefix}.${segid}.tmp.pdb"
        lappend seg_pdbs $seg_pdb
        [atomselect top "segid $segid"] writepdb $seg_pdb

        segment $segid {
            pdb $seg_pdb
            if {$segid == $mutate_segid} {
                mutate $mutate_resid $mutate_to
            }
         }
         coordpdb $seg_pdb $segid
    }
    guesscoord
    writepsf "${out_prefix}.psf"
    writepdb "${out_prefix}.pdb"

    # Clean up temporary files
    foreach seg_pdb $seg_pdbs {
        file delete $seg_pdb
    }
}

proc _do_4term_extrabonds_line {molid params_dict keyword {fd stdout}} {
    set all [atomselect top all]
    foreach lipid [dict keys $params_dict] {
        foreach segid [lunique [$all get segid]] {
            set sel [atomselect $molid "segid $segid and resname $lipid"]
            foreach resid [lunique [$sel get resid]] {
                # Note that lrange is inclusive in its range
                set atomnames [lrange [dict get $params_dict $lipid] 0 3]
                set ref_val [lindex [dict get $params_dict $lipid] 4]
                set res_sel [atomselect $molid "segid $segid and resname $lipid and resid $resid and name $atomnames"]
                set indexes [$res_sel get index]
                puts $fd "$keyword $indexes \$FC $ref_val"

                $res_sel delete
            }
            $sel delete
        }
    }
    $all delete
}

# Generates lipid restraints as used by the CHARMM-GUI lipid equilibration protocol.
# Useful when you are, for example, doing AFEP on a lipid that CHARMM-GUI generated.
proc generate_lipid_restraints {molid} {
    # First we deal with the head atoms, and generate ${lipid}_head_{upper,lower}.ref files
    # that specify what atoms are to be restrained in the lipid heads.

    # Tcl gives a useless error message if you have a space after the line continuation backslash
    set head_atoms_by_lipid [dict create \
        "POEA" {P N "C1.*" "O1.*"} \
        "POPE" {P N "C1.*" "O1.*"} \
        "POPA" {P C1 "O1.*"} \
    ]

    puts "generate_lipid_restraints: Warning: I only know heads of:\n    [lsort [dict keys $head_atoms_by_lipid]]\nYou _must_ review the source because this proc is shady!"
    set all [atomselect $molid all]

    # Mark head atoms in upper (z > 0) and lower (z < 0) leaflets
    foreach lipid [dict keys $head_atoms_by_lipid] {
        $all set beta 0
        set head_atoms [dict get $head_atoms_by_lipid $lipid]
        set sel [atomselect $molid "resname $lipid and z > 0 and name $head_atoms"]
        $sel set beta 1
        $all writepdb "[string tolower $lipid]_head_upper.ref"
        $sel delete

        $all set beta 0
        set sel [atomselect $molid "resname $lipid and z < 0 and name $head_atoms"]
        $sel set beta 1
        $all writepdb "[string tolower $lipid]_head_lower.ref"
        $sel delete
    }

    # Now generate dihe.txt, which is used as input to NAMD extraBonds, and consists of lines like
    # DIHEDRAL a b c d $FC 0.0
    # IMPROPER a b c d $FC 120.0
    # where a b c d are atom indices (zero-based!) in the lipid tails.
    # This will be different according to lipid type, and are encoded in an inscrutable CHARMM script.
    # So I'm only including certain lipid types piecemeal.
    # First four are atom names, last is reference value
    set eb_DIHEDRAL [dict create \
        "POPE" "C28 C29 C210 C211 0.0" \
    ]
    set eb_IMPROPER [dict create \
        "POPE" "C3 C1 C2 O21 120.0" \
    ]

    set fd [open "dihe.txt" "w"]
    _do_4term_extrabonds_line $molid $eb_DIHEDRAL "DIHEDRAL" $fd
    _do_4term_extrabonds_line $molid $eb_IMPROPER "IMPROPER" $fd
    close $fd
    $all delete
}


# Generates (rough) table of corresponding residues between structurally similar structures.
# Useful for e.g. WSMuOR vs MuOR - which residues in the derived version correspond to the original?
#   distance_cutoff: how far away is too far between CA atoms to say the residues correspond
proc guess_corresponding_residues {molid segid ref_molid {distance_cutoff 5}} {
    # The structures must first be aligned
    align_trajectory $molid $ref_molid

    # Iterate through each residue of the molecule whose residues we are curious about
    # Since each residue has an alpha carbon, use that.
    # TODO: Does the segid really matter? Can we just get rid of it?
    set mol_ca_sel [atomselect $molid "segid $segid and protein and name CA"]

    set atom_indices [$mol_ca_sel list]
    # Since atom indices and resids are simply labels and nothing can be assumed
    # about their ordering or contiguity, we use a separate index to iterate through
    # both simultaneously
    foreach atom_idx $atom_indices {
        # set atom_idx [lindex $atom_indices $i]
        # What's the nearest atom in the reference molecule to this one?
        # First, we need the coordinates of this atom
        set this_atom_sel [atomselect $molid "index $atom_idx"]
        set this_resid [$this_atom_sel get resid]
        set x [$this_atom_sel get x]
        set y [$this_atom_sel get y]
        set z [$this_atom_sel get z]
        set xyz [lindex [$this_atom_sel get {x y z}] 0]
        # Find all CA atoms within $distance_cutoff of atom atom_idx
        # Conceivably, we could get all the coordinates of everything and do the N^2 distance matrix calculation.
        # However, the thought of doing anything in Tcl involving some sort of data structure is less than fun.
        set d2 [expr $distance_cutoff * $distance_cutoff]
        set nearest_ref_ca_sel [atomselect $ref_molid "protein and name CA and ((x-$x)*(x-$x) + (y-$y)*(y-$y) + (z-$z)*(z-$z)) < $d2"]
        # We have an atomselect of the nearby alpha carbons in the reference molecule, but it isn't ordered.
        # So we will need to calculate all the distances.
        set nearest_ref_ca_idx [$nearest_ref_ca_sel get index]
        set nearest_ref_ca_xyz [$nearest_ref_ca_sel get {x y z}]
        # Now iterate through each and find the smallest distance
        set num_nearest [llength $nearest_ref_ca_idx]
        # Couldn't find anything within cutoff? Then this residue is unsalvageable
        if {$num_nearest <= 0} {
            puts "guess_corresponding_residues: No corresponding residue in reference for resid $this_resid"
            continue
        }
        set closest_j -1
        set closest_dist [expr $distance_cutoff + 1000]
        for {set j 0} {$j < $num_nearest} {incr j} {
            # Calculate distance and save the relevant atom if it was the closest
            set dist [vecdist $xyz [lindex $nearest_ref_ca_xyz $j]]
            if {$dist < $closest_dist} {
                set closest_j $j
                set closest_dist $dist
            }
        }
        set closest_resid [lindex [$nearest_ref_ca_sel get resid] $closest_j]
        set closest_resname [lindex [$nearest_ref_ca_sel get resname] $closest_j]
        set this_resname [$this_atom_sel get resname]
        puts "$this_resname$this_resid -> $closest_resname$closest_resid ($closest_dist A)"

        # Garbage collection is for wimps and better interpreted languages
        $nearest_ref_ca_sel delete
        $this_atom_sel delete
    }
    $mol_ca_sel delete
}
