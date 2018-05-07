# Calculates RMSD only of binding pocket residues between water-soluble mu opioid reeptor
#
# Example:
#   muor_binding_pocket_rmsd 3 0 -42
# To use molid 3 as ref, molid 0 as WSMuOR_TMv2 trajectory, and -42 as offset
# For WSMuOR_TMv1, offset should be -62
# set REF_BINDING_POCKET_RESIDUES {233 148 297 300 296 151 326 293 147}

proc muor_binding_pocket_rmsd {ref_molid molid {residue_offset -42} {ref_residues {233 148 297 300 296 151 326 293 147}}} {
	# Add offset to resids in the real MuOR to get the corresponding offset
	# in the WSMuOR
	set binding_pocket_residues {}
	set ref_res_sel {}
	set ws_res_sel {}
	set column_labels {}
	
	# Cache atomselect objects for each residue
	foreach resid $ref_residues {
		set ws_resid [expr $resid + $residue_offset]
		lappend binding_pocket_residues $ws_resid
		lappend ws_res_sel [atomselect $molid "backbone and resid $ws_resid"]
		# Am I living in a dream world? I have to use a command called "lindex" instead of a sane
		# array subscripting notation, like [] as in essentially every other language?
		set ref_res [atomselect $ref_molid "backbone and resid $resid"]
		set ref_res_name [string totitle [lindex [$ref_res get resname] 0]]
		lappend ref_res_sel $ref_res
		lappend column_labels "$ref_res_name$resid"
		
	}

	# Dump this into a file
	set fd [open "binding_pocket_rmsd.csv" "w"]

	# Print CSV column headers
	lappend column_labels "Combined"
	puts $fd [join $column_labels ","]

	# Select the atoms we're going to use.
	# We do this outside the loop because each time you do this, a new object is created
	set ref [atomselect $ref_molid "backbone and resid $ref_residues"]
	set ref_all [atomselect $ref_molid all]
	set ws [atomselect $molid "backbone and resid $binding_pocket_residues"]
	set ws_all [atomselect $molid all]
	set num_residues [llength $binding_pocket_residues]

	
	# Iterate over frames
	set num_frames [molinfo $molid get numframes]
	for {set frame 0} {$frame < $num_frames} {incr frame} {
		# Set the current frame
		$ws frame $frame
		# Do the alignment
		set trans_mat [measure fit $ws $ref]
		$ws move $trans_mat
		set overall_rmsd [measure rmsd $ws $ref]
		set per_res_rmsd {}
		for {set i 0} {$i < $num_residues} {incr i} {
			[lindex $ws_res_sel $i] frame $frame
			lappend per_res_rmsd [format "%.2f" [measure rmsd [lindex $ref_res_sel $i] [lindex $ws_res_sel $i]]]
		}
		# Print the individual residue RMSDs
		puts -nonewline $fd [join $per_res_rmsd ","]
		# Print the overall RMSD of all residues in question
		puts -nonewline $fd ","
		puts $fd "$overall_rmsd"
	}

	# Close the output file
	close $fd

	# Delete all those atomselect objects because apparently there is no garbage collection
	$ref delete
	$ws delete
	foreach sel $ref_res_sel { $sel delete }
	foreach sel $ws_res_sel { $sel delete }
}
