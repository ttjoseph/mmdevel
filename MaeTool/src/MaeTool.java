import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

public class MaeTool
{
	public static void main(String[] args)
	{
		System.err.println("MaeTool - by Tom Joseph <thomas.joseph@mssm.edu>\n");
		//main_makeFullSystem(args);
		//main_deleteSomeAtoms(args);
		//main_addBSC0Dihedrals(args);
		//main_fixAtomicNumbersForVMD(args);
		//main_testFEPSetup(args);
		main_convertToPDB(args);
	}
	
	public static void main_testFEPSetup(String[] args)
	{
		final String baseDir = "/compass/users1/tom/PertTest/";
		MaeNode orig = null, pert = null;
		try
		{
			orig = MaeNode.fromCMSFile(baseDir + "original.cms");
			pert = MaeNode.fromCMSFile(baseDir + "perturbed.cms");
		} catch (IOException e)
		{
			e.printStackTrace();
			return;
		}
		
		try
		{
			//orig.setupForAlchemicalFEP(pert);
			//orig.toStream(System.out);
			orig.getNode("f_m_ct").getNode("m_atom").dataRecordsAsSQLToStream(System.out);
			throw new MaeNodeException("Shut up");
		} catch(MaeNodeException e)
		{
			e.printStackTrace();
			return;
		}		
	}
	
	public static void main_fixAtomicNumbersForVMD(String[] args)
	{
		final String baseDir = "/compass/users1/tom/TtAgoAA6/";
		MaeNode cms = null;
		try
		{
			cms = MaeNode.fromCMSFile(baseDir + "TtAA6_1-in.cms");
		} catch (IOException e)
		{
			e.printStackTrace();
			return;
		}
		
		try
		{
			fixAtomicNumbersForVMD(cms);
			cms.toStream(System.out);
		} catch(MaeNodeException e)
		{
			e.printStackTrace();
			return;
		}
	}
	
	public static void main_deleteSomeAtoms(String[] args)
	{
		final String baseDir = "/compass/users1/tom/TtAgoAA6/Scratch/Rebuild/";
		MaeNode cms = null;
		try
		{
			cms = MaeNode.fromCMSFile(baseDir + "soluteAA6.mae");
		} catch (IOException e)
		{
			e.printStackTrace();
			return;
		}
		
		try
		{
			MaeNode ct = cms.getNode("f_m_ct");
			//ct.shiftResidueIDs(700, "MD6");
			//ct.shiftResidueIDs(800, "DA", "DG", "DT", "RA", "RC", "RG", "RU");
			//ct.deleteAtom(11372);
			//ct.deleteAtom(11371);
			//ct.deleteAtom(11370);
			cms.toStream(System.out);
			throw new MaeNodeException("Shut up, compiler");
		} catch(MaeNodeException e)
		{
			e.printStackTrace();
			return;
		}		
	}
	
	public static void main_convertToPDB(String[] args)
	{
		System.err.println("Converting your mae/cms file to PDB.\n" +
				"Remember, you may still need to massage it to get other programs to accept it.");
		
		final String baseDir = "/compass/users1/tom/TtAgoDRtoAA6/Build/";
		MaeNode cms = null;
		try
		{
			cms = MaeNode.fromCMSFile(baseDir + "TtAgoDR.cms");
		} catch (IOException e)
		{
			e.printStackTrace();
			return;
		}
		
		//MaeNode m_atom = cms.getNode("f_m_ct").getNode("m_atom");
		//writeMinimalPSF(System.out, m_atom, "CA");
		try
		{
			cms.soluteOnlyToStreamAsPDB(System.out);
		} catch (MaeNodeException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public static void main_addBSC0Dihedrals(String[] args)
	{
		final String baseDir = "/compass/users1/tom/TtAgoAA6/Scratch/Rebuild/";
		MaeNode cms = null;
		try
		{
			cms = MaeNode.fromCMSFile(baseDir + "allAA6-almost.cms");
		} catch (IOException e)
		{
			e.printStackTrace();
			return;
		}
		
		try
		{
			addBSC0Dihedrals(cms);
			cms.toStream(System.out);
		} catch(MaeNodeException e)
		{
			e.printStackTrace();
			return;
		}
		
	}
	
	public static void main_makeFullSystem(String[] args)
	{
		MaeNode solute = null, solvent = null, ions = null;
		final String baseDir = "/compass/users1/tom/TtAgoAA6/Scratch/Rebuild/";
		try
		{
			solute = MaeNode.fromCMSFile(baseDir + "soluteAA6-almost.cms");
			solvent = MaeNode.fromCMSFile(baseDir + "solvent.mae");
			ions = MaeNode.fromCMSFile(baseDir + "na.mae");
		} catch (IOException e)
		{
			e.printStackTrace();
			return;
		}

		try
		{
			// TODO: FF terms are in solute, not fullsystem?
			// TODO: BSC0 stuff doesn't work on OPLS-paramterized system, must run through Viparr first
			MaeNode fullSystem = makeFullSystemCT(
					solute.getNode("f_m_ct"), 
					solvent.getNode("f_m_ct"),
					ions.getNode("f_m_ct"));
			
			// Assemble a complete CMS file from these pieces
			
			// Set s_ffio_ct_type and stuff correctly
			solute.getNode("f_m_ct").setString("s_ffio_ct_type", "solute");
			ions.getNode("f_m_ct").setString("s_ffio_ct_type", "ion");

			// We'll assume that the solvent has the PBC parameters
			solute.getNode("f_m_ct").copyPBCParameters(solvent.getNode("f_m_ct"));
			ions.getNode("f_m_ct").copyPBCParameters(solvent.getNode("f_m_ct"));
						
			String[] nothing = {};
			MaeNode all = new MaeNode("HEAD", false, nothing);
			// Gotta have the header block, so we copy an existing one
			all.addChildNode(solute.getNode("anonymous"));
			all.addChildNode(fullSystem);
			all.addChildNode(solute.getNode("f_m_ct"));
			all.addChildNode(solvent.getNode("f_m_ct"));
			all.addChildNode(ions.getNode("f_m_ct"));
			
			all.toStream(System.out);
		} catch (MaeNodeException e)
		{ e.printStackTrace(); }
	}

	public static void main_convertToPSF(String[] args)
	{
		final String baseDir = "/compass/users1/tom/TtAgoAA6/";
		MaeNode cms = null;
		try
		{
			cms = MaeNode.fromCMSFile(baseDir + "TtAgoAA6.fixed.cms");
		} catch (IOException e)
		{
			e.printStackTrace();
			return;
		}
		
		//MaeNode m_atom = cms.getNode("f_m_ct").getNode("m_atom");
		//writeMinimalPSF(System.out, m_atom, "CA");
		try
		{
			cms.soluteOnlyToStreamAsPSF(System.out);
		} catch (MaeNodeException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * VMD apparently doesn't like atomic numbers >= 200, so this changes them to 1.
	 * 
	 * @param cms
	 * @throws MaeNodeException
	 */
	public static void fixAtomicNumbersForVMD(MaeNode cms) throws MaeNodeException
	{
		System.err.println("Fixing MD6 atomic numbers so VMD won't barf.");
		
		// VMD chokes on nonstandard atomic numbers for some reason, so change
		// MD6 atomic numbers
		// from 200-203 to 1.
		// The error message talks about namelists and has nothing to do with
		// the actual
		// reason for failure.
		for (int i = 0; i < 4; i++)
		{
			MaeNode m_atom = cms.getNode("f_m_ct", i).getNode("m_atom");
			List<String> atomicNumbers = m_atom.getColumn("i_m_atomic_number");
			for (int j = 0; j < atomicNumbers.size(); j++)
			{
				String num = atomicNumbers.get(j);
				if (Integer.parseInt(num) >= 200)
					num = "1";
				atomicNumbers.set(j, num);
			}
		}
	}

	public static MaeNode addBSC0Dihedrals(MaeNode cms) throws MaeNodeException
	{
		// System.err.println("Format version: " +
		// cms.getNode("anonymous").getString("s_m_m2io_version", 0));

		ArrayList<ProperDihedralTerm> terms = new ArrayList<ProperDihedralTerm>();

		// ParmBSC0 terms, AMBER-style
		terms.add(new ProperDihedralTerm("OS", "P", "OS", "CI", 1, 1.256531, 351.95960, 2));
		terms.add(new ProperDihedralTerm("OS", "P", "OS", "CI", 1, 0.354858, 357.24748, 3));
		terms.add(new ProperDihedralTerm("OH", "P", "OS", "CI", 1, 1.256531, 351.95960, 2));
		terms.add(new ProperDihedralTerm("OH", "P", "OS", "CI", 1, 0.354858, 357.24748, 3));
		terms.add(new ProperDihedralTerm("CT", "CT", "CI", "OS", 1, 0.092102, 295.63279, 2));
		terms.add(new ProperDihedralTerm("CT", "CT", "CI", "OS", 1, 0.962830, 348.09535, 3));
		terms.add(new ProperDihedralTerm("CT", "CT", "CI", "OH", 1, 0.092102, 295.63279, 2));
		terms.add(new ProperDihedralTerm("CT", "CT", "CI", "OH", 1, 0.962830, 348.09535, 3));

		MaeNode soluteCT = cms.getCTNodeByType("solute");
		MaeNode ffio_sites = soluteCT.getNode("ffio_ff").getNode("ffio_sites");
		MaeNode ffio_dihedrals = soluteCT.getNode("ffio_ff").getNode("ffio_dihedrals");

		for (ProperDihedralTerm term : terms)
			applyDihedralTerm(ffio_sites, ffio_dihedrals, term);

		return cms;
	}

	/**
	 * Applies a single proper dihedral term to the force field nodes of a CT block.
	 * 
	 * @param ffio_sites
	 * @param ffio_dihedrals
	 * @param dt
	 */
	public static void applyDihedralTerm(MaeNode ffio_sites,
			MaeNode ffio_dihedrals, ProperDihedralTerm dt)
	{
		System.err.println("Applying dihedral term: " + dt.toString());
		Set<ProperDihedralTerm.FourAtoms> atomIDs = dt.getAtomIndices(ffio_sites, ffio_dihedrals);
		dt.applyToCMS(ffio_dihedrals, atomIDs);
	}
	
	/**
	 * Concatenates the m_atom and m_bond child blocks of CT blocks given,
	 * adjusting the m_bond child blocks to reflect the changed atom indices.
	 * This yields a full_system CT block that you probably want to copy and
	 * paste into a CMS file.
	 * 
	 * @param ct_blocks
	 * @return
	 */
	public static MaeNode makeFullSystemCT(MaeNode... ct_blocks)
			throws MaeNodeException
	{
		String[] headings = { "s_m_title", "r_chorus_box_ax", "r_chorus_box_ay",
				"r_chorus_box_az", "r_chorus_box_bx", "r_chorus_box_by","r_chorus_box_bz",
				"r_chorus_box_cx", "r_chorus_box_cy", "r_chorus_box_cz", "s_ffio_ct_type" };

		MaeNode fullSystem = new MaeNode("f_m_ct", false, headings);

		// We need to merge the m_atom and m_bond blocks from each ct block.
		// This entails renumbering the atoms, but not the residue numbers
		// The m_bond entries must be updated with the new atom numbers of its
		// corresponding m_atom block.

		HashMap<String, String> atomRecord = new HashMap<String, String>();
		// ArrayList<String> atomHeadings = new ArrayList<String>();
		String[] atomHeadings = { "i_m_mmod_type", "r_m_x_coord",
				"r_m_y_coord", "r_m_z_coord", "i_m_residue_number",
				"s_m_insertion_code", "s_m_mmod_res", "s_m_chain_name",
				"i_m_color", "r_m_charge1", "r_m_charge2",
				"s_m_pdb_residue_name", "s_m_pdb_atom_name", "s_m_grow_name",
				"i_m_atomic_number", "i_m_formal_charge", "s_m_atom_name",
				"i_m_secondary_structure", "i_m_Hcount",
				"i_m_pdb_convert_problem", "i_ppw_het" };
		atomRecord.put("s_m_label_format", "1"); // atomHeadings.add("s_m_label_format");
		atomRecord.put("i_m_label_color", "1"); // atomHeadings.add("i_m_label_color");
		// atomHeadings.add("s_m_label_user_text");
		atomRecord.put("r_m_pdb_tfactor", "1"); // atomHeadings.add("r_m_pdb_tfactor");
		atomRecord.put("r_m_pdb_occupancy", "1"); // atomHeadings.add("r_m_pdb_occupancy");
		atomRecord.put("i_m_representation", "0"); // atomHeadings.add("i_m_representation");
		atomRecord.put("i_m_visibility", "1"); // atomHeadings.add("i_m_visibility");
		atomRecord.put("i_m_template_index", "0"); // atomHeadings.add("i_m_template_index");

		HashMap<String, String> bondRecord = new HashMap<String, String>();
		String[] bondHeadings = { "i_m_from", "i_m_to", "i_m_order" };

		MaeNode newAtoms = new MaeNode("m_atom", true, atomHeadings);
		MaeNode newBonds = new MaeNode("m_bond", true, bondHeadings);

		// Process CT blocks, extracting m_atom and m_bond blocks
		int atomIndexOffset = 0;
		for (MaeNode ct : ct_blocks)
		{
			// Does this node contain PBC parameters? If so, copy them.
			// XXX: This means you will end up with whatever the last PBCs were.
			// Normally they are the same in every CT block for a given system but
			// who knows what the future may hold...
			fullSystem.copyPBCParameters(ct);
			
			MaeNode m_atom = ct.getNode("m_atom");
			MaeNode m_bond = ct.getNode("m_bond");
			if (m_atom.isDataArray() == false)
			{
				System.err.println("All m_atom blocks should be data arrays.  I hope you know what that means.");
				System.exit(0);
			}

			//for (String heading : m_atom.getHeadings())
			//	System.err.println("makeFullSystemCT: m_atom has heading " + heading);

			for (int i = 0; i < m_atom.numDataRecords(); i++)
			{
				for (String heading : atomHeadings)
				{
					String value = m_atom.getString(heading, i);
					if (value == null) // Not all m_atom blocks have the same set of headings
						value = "<>"; // I think this is the "no data" token for
										// mae files
				
					atomRecord.put(heading, value);
				}
				
				// Put solvent molecules in a different chain to placate Viparr
				// Just sort of picking these chain names at random
				String pdbName = atomRecord.get("s_m_pdb_residue_name").trim();

				if(pdbName.equals("T3P"))
					atomRecord.put("s_m_chain_name", "S");
				/*
				else if(pdbName.equals("MD6"))
				{
					System.err.println("Adding 700 to MD6 residue number.");
					// Offset residue numbers so Viparr doesn't whine about it
					atomRecord.put("i_m_residue_number",
							Integer.toString(Integer.parseInt(atomRecord.get("i_m_residue_number")) + 700));
				}
				else if(pdbName.equals("DA")
						|| pdbName.equals("DG")
						|| pdbName.equals("DT")
						|| pdbName.equals("RA")
						|| pdbName.equals("RC")
						|| pdbName.equals("RG")
						|| pdbName.equals("RU"))
				{
					// Offset residue numbers so Viparr doesn't whine about it
					atomRecord.put("i_m_residue_number",
							Integer.toString(Integer.parseInt(atomRecord.get("i_m_residue_number")) + 800));
				}
				*/
				// I guess we don't have to massage this data like we do for
				// m_bond
				
				// copy from old m_atom to new m_atom
				newAtoms.addDataRecord(atomRecord);
				
				// Update old m_atom with whatever changes we made
				m_atom.setDataRecord(atomRecord, i);
			}
			
			System.err.println("makeFullSystemCT: T3P residues are marked as chain S.  Hope that's OK.");

			// m_bond may not exist, like in CT blocks composed exclusively of
			// point-charge ions
			if (m_bond != null)
			{
				for (int i = 0; i < m_bond.numDataRecords(); i++)
				{
					for (String heading : bondHeadings)
					{
						String value = m_bond.getString(heading, i);
						if (value == null)
							value = "<>"; // See above
						bondRecord.put(heading, value);
					}

					// Update atom references so they are consistent with the
					// new m_atom block
					int fromAtom = Integer.parseInt(bondRecord.get("i_m_from"));
					int toAtom = Integer.parseInt(bondRecord.get("i_m_to"));
					fromAtom += atomIndexOffset;
					toAtom += atomIndexOffset;
					bondRecord.put("i_m_from", Integer.toString(fromAtom));
					bondRecord.put("i_m_to", Integer.toString(toAtom));
					newBonds.addDataRecord(bondRecord);
				}
			}

			atomIndexOffset += m_atom.numDataRecords();
		}

		fullSystem.addChildNode(newAtoms);
		fullSystem.addChildNode(newBonds);
		fullSystem.setString("s_m_title", "full system");
		fullSystem.setString("s_ffio_ct_type", "full_system");

		//System.err.println("makeFullSystemCT: Remember, you still have to manually add PBC parameters!");
		return fullSystem;
	}

	/**
	 * Writes a PSF file that only contains atom indices of atoms that have a
	 * specific name. This is useful for feeding to postprocessing programs, not
	 * for actual simulation. For example, essential dynamics. You can use VMD
	 * to make a corresponding DCD trajectory.
	 * 
	 * @param out
	 * @param m_atom
	 * @param atomName
	 */
	public static void writeMinimalPSF(java.io.PrintStream out, MaeNode m_atom,
			String queryAtomName)
	{
		// Is this really an atom block?
		if (m_atom.getName().equals("m_atom") == false
				|| m_atom.isDataArray() == false)
		{
			System.err
					.println("writeMinimalPSF: This MaeNode doesn't look like an m_atom block.  I give up.");
			return;
		}

		// We'll just examine the PDB name column
		ArrayList<String> atomNames = m_atom.getColumn("s_m_pdb_atom_name");
		if (atomNames == null)
		{
			System.err
					.println("writeMinimalPSF: This m_atom block doesn't seem to have the s_m_pdb_atom_name field.  I need that.");
			return;
		}

		// Loop through the atom block, find atoms with names matching atomName,
		// and add their
		// indices to a list
		ArrayList<Integer> atomIDs = new ArrayList<Integer>();
		for (int i = 0; i < m_atom.numDataRecords(); i++)
		{
			if (atomNames.get(i).trim().equals(queryAtomName))
				atomIDs.add(i);
		}

		// PSF global header
		out.println("PSF CMAP\n");
		out.println("       2 !NTITLE");
		out.println(" REMARKS Generated by a program written by Tom Joseph");
		out.println(" REMARKS <thomas.joseph@mssm.edu>\n");

		// Atom block
		out.printf("%1$8d !NATOM\n", atomIDs.size());
		// 1231 !NATOM
		// 1 U 1 MET N NH3 -0.300000 14.0070 0

		int printedAtomID = 1;
		for (int atomID : atomIDs)
		{
			String segID = m_atom.getString("s_m_chain_name", atomID).trim();
			String resID = m_atom.getString("i_m_residue_number", atomID)
					.trim();
			String resName = m_atom.getString("s_m_pdb_residue_name", atomID)
					.trim();
			String atomName = m_atom.getString("s_m_pdb_atom_name", atomID)
					.trim();

			// dump the atom info
			// TODO: what is the correct format here?
			// atom ID, seg ID, res ID, res name, atom name, type, charge, mass,
			// 0
			// We don't really care about type, charge, mass
			out
					.printf(
							"%1$8d %2$-4s %3$-4d %4$-4s %5$-4s %6$-4s %7$12.6f %8$10.4f\n",
							printedAtomID, "X", Integer.parseInt(resID),
							resName, atomName, atomName, 0.0, 1.0);

			printedAtomID++;
		}

		out.println();
		// print dummy placeholder blocks for bonds, angles, dihedrals, etc
		out.printf("%1$8d !NBOND\n\n", 0);
		out.printf("%1$8d !NTHETA\n\n", 0);
		out.printf("%1$8d !NPHI\n\n", 0);
		out.printf("%1$8d !NIMPHI\n\n", 0);
		out.printf("%1$8d !NCRTERM\n\n", 0);

	}
}
