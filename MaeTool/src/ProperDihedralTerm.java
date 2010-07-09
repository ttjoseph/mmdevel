import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class ProperDihedralTerm
{
	/**
	 * Encapsulates a set of four atom indices meant for a dihedral term. They
	 * can be in forward or reverse order, and this class is aware of that and
	 * the behavior of the comparator and such reflect this.
	 */
	public class FourAtoms implements Comparator<FourAtoms>
	{
		public int i, j, k, l;

		public FourAtoms(int i_, int j_, int k_, int l_)
		{
			i = i_;
			j = j_;
			k = k_;
			l = l_;
		}

		public int compare(FourAtoms o1, FourAtoms o2)
		{
			int ai = o1.i, aj = o1.j, ak = o1.k, al = o1.l;
			int bi = o2.i, bj = o2.j, bk = o2.k, bl = o1.l;

			// Ensure that i <= l for both objects, because the order
			// of the elements can be reversed and the parent objects
			// would still be considered equal, as a consequence of
			// applying to a dihedral term. Doing this removes ambiguity.
			if (ai > al)
			{
				ai = o1.l;
				aj = o1.k;
				ak = o1.j;
				bl = o1.i;
			}
			if (bi > bl)
			{
				bi = o2.l;
				bj = o2.k;
				bk = o2.j;
				bl = o2.i;
			}

			// Do the actual comparison
			if (ai < bi)
				return -1;
			if (aj < bj)
				return -1;
			if (ak < bk)
				return -1;
			if (al < bl)
				return -1;

			// Since o1 is definitely not "less than" o2, it might be equal
			if (ai == bi && aj == bj && ak == bk && al == bl)
				return 0;

			// Not less than or equal, so o2 must be greater
			return 1;
		}
	}

	/** Names of the atoms that this dihedral term applies to. */
	private String atom1_, atom2_, atom3_, atom4_;
	/** Proper dihedral parameters, in Viparr Proper_Trig format. */
	private double phase_, fc_[];

	/**
	 * Constructs dihedral term from AMBER DIHE-style parameters (as in frcmod
	 * files).
	 * 
	 * @param atom1
	 * @param atom2
	 * @param atom3
	 * @param atom4
	 * @param idivf
	 * @param pk
	 * @param phase
	 *            Phase, in degrees
	 * @param pn
	 */
	public ProperDihedralTerm(String atom1, String atom2, String atom3,
			String atom4, double idivf, double pk, double phase, int pn)
	{
		atom1_ = atom1;
		atom2_ = atom2;
		atom3_ = atom3;
		atom4_ = atom4;
		phase_ = phase;
		fc_ = new double[7];
		fc_[pn] = pk / idivf;

		for (double fc : fc_)
			fc_[0] += Math.abs(fc);
	}

	/**
	 * Checks whether this dihedral term applies to the four atom type names
	 * specified.
	 * 
	 * @param atom1
	 * @param atom2
	 * @param atom3
	 * @param atom4
	 * @return True if it does, false if not.
	 */
	public boolean appliesTo(String atom1, String atom2, String atom3,
			String atom4)
	{
		// Can be in forward or reverse order
		return ((atom1.equals(atom1_) && atom2.equals(atom2_)
				&& atom3.equals(atom3_) && atom4.equals(atom4_)) || (atom4
				.equals(atom1_)
				&& atom3.equals(atom2_) && atom2.equals(atom3_) && atom1
				.equals(atom4_)));
	}

	public String toString()
	{
		StringBuffer s = new StringBuffer("\"" + atom1_ + "\" \"" + atom2_
				+ "\" \"" + atom3_ + "\" \"" + atom4_ + "\" " + phase_);
		for (double fc : fc_)
			s.append(" " + fc);
		return s.toString();
	}

	/**
	 * Returns a set of atom indices that this proper dihedral term applies to.
	 * Presumably you will take this set of atom indices and pass them on to
	 * applyToCMS. Yes, this is a silly way to design it.
	 * 
	 * @param m_atom
	 *            The m_atom block from the mae file
	 * @param ffio_dihedrals
	 *            The ffio_dihedrals block from the mae file
	 */
	public Set<FourAtoms> getAtomIndices(MaeNode ffio_sites,
			MaeNode ffio_dihedrals)
	{
		// This dihedral applies to a set of four atom names, in either forward
		// or reverse order.
		// We don't know which atom indices these correspond to. We need to find
		// valid sets of
		// atom indices that reflect atoms that are actually bonded in the
		// structure. Happily,
		// we can just look for existing dihedral terms that share the same atom
		// names as us.

		HashSet<FourAtoms> set = new HashSet<FourAtoms>();

		int num_dihedrals = ffio_dihedrals.numDataRecords();
		System.err.println("getAtomIndices: Examining " + num_dihedrals
				+ " dihedral terms");
		for (int i = 0; i < num_dihedrals; i++)
		{
			int ai = ffio_dihedrals.getInteger("i_ffio_ai", i);
			int aj = ffio_dihedrals.getInteger("i_ffio_aj", i);
			int ak = ffio_dihedrals.getInteger("i_ffio_ak", i);
			int al = ffio_dihedrals.getInteger("i_ffio_al", i);
			String type = ffio_dihedrals.getString("s_ffio_funct", i);
			// Atom IDs are 1-based in the file but our API expects them to be
			// 0-based
			String atom1 = ffio_sites.getString("s_ffio_vdwtype", ai - 1)
					.trim();
			String atom2 = ffio_sites.getString("s_ffio_vdwtype", aj - 1)
					.trim();
			String atom3 = ffio_sites.getString("s_ffio_vdwtype", ak - 1)
					.trim();
			String atom4 = ffio_sites.getString("s_ffio_vdwtype", al - 1)
					.trim();

			// if(i == 118)
			// {
			// System.err.println("getAtomIndices: " + ai + " " + aj + " " + ak
			// + " " + al);
			// System.err.println("getAtomIndices: " + atom1 + " " + atom2 + " "
			// + atom3 + " " + atom4);
			// }
			if (type.equals("Proper_Trig")
					&& this.appliesTo(atom1, atom2, atom3, atom4))
			{
				System.err.println("YES! Dihedral " + (i + 1)
						+ " applies to the same atoms I do!");
				FourAtoms atomIDs = new FourAtoms(ai, aj, ak, al);
				if (set.contains(atomIDs) == false)
					set.add(atomIDs);
			}
		}

		return set;
	}

	/**
	 * Applies this dihedral to a set of atom indices, which you probably got
	 * from getAtomIndices. You can then output this CMS file and hopefully use
	 * it with Desmond.
	 * 
	 * @param atomIDs
	 *            Sets of atom IDs that this dihedral term applies to.
	 * @param ffio_dihedrals
	 */
	public void applyToCMS(MaeNode ffio_dihedrals, Set<FourAtoms> atomIDs)
	{
		// For each set of 4 atoms, add this dihedral term using those IDs to
		// the
		// specified ffio_dihedrals block, which presumably is part of a CMS
		HashMap<String, String> record = new HashMap<String, String>();
		record.put("s_ffio_funct", "Proper_Trig");
		record.put("r_ffio_c0", Double.toString(phase_)); // c0 is actually the
															// phase
		record.put("r_ffio_c1", Double.toString(fc_[0]));
		record.put("r_ffio_c2", Double.toString(fc_[1]));
		record.put("r_ffio_c3", Double.toString(fc_[2]));
		record.put("r_ffio_c4", Double.toString(fc_[3]));
		record.put("r_ffio_c5", Double.toString(fc_[4]));
		record.put("r_ffio_c6", Double.toString(fc_[5]));
		record.put("r_ffio_c7", Double.toString(fc_[6]));

		for (FourAtoms atoms : atomIDs)
		{
			// Actually add this dihedral...
			record.put("i_ffio_ai", Integer.toString(atoms.i));
			record.put("i_ffio_aj", Integer.toString(atoms.j));
			record.put("i_ffio_ak", Integer.toString(atoms.k));
			record.put("i_ffio_al", Integer.toString(atoms.l));
			
			ffio_dihedrals.addDataRecord(record);
		}

	}
}
