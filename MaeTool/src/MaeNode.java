import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.StreamTokenizer;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.ArrayList;

/**
 * Encapsulates the mae/maeff/cms file format used in Maestro and Desmond. A
 * MaeNode may have child MaeNodes.
 * 
 * @author Tom Joseph <thomas.joseph@mssm.edu>
 */
public class MaeNode
{
	/** List of child nodes of this node */
	private TreeMap<String, ArrayList<MaeNode>> children_ = new TreeMap<String, ArrayList<MaeNode>>();
	/**
	 * List of heading names - that is, the keys of children_. This variable is
	 * here because the heading names have an ordering that must be enforced.
	 */
	private ArrayList<String> headings_ = new ArrayList<String>();

	/** Represents a column of the table contained in this node */
	private TreeMap<String, ArrayList<String>> dataColumns_ = new TreeMap<String, ArrayList<String>>();
	private int numDataRecords_ = 0;

	/** A node without a name is given the name "anonymous" */
	private String name_ = "anonymous";
	private boolean isDataArray_ = false;

	/** Indent level, for making printing to a stream pretty */
	private static String indent_ = "";

	/**
	 * Increases the indent level for printing. Purely of cosmetic importance.
	 */
	private static void increaseIndent()
	{
		indent_ = indent_ + "  ";
	}

	/**
	 * Decreases the indent level for printing. Purely cosmetic.
	 */
	private static void decreaseIndent()
	{
		if (indent_.length() >= 2)
			indent_ = indent_.substring(2);
	}

	/**
	 * Returns the node of the given name. Since there can be more than one node
	 * with the same name, you also have to specify an index. The first index is
	 * 0.
	 * 
	 * @param name
	 * @param index
	 * @return The node you want, or null if not found
	 */
	public MaeNode getNode(String name, int index)
	{
		ArrayList<MaeNode> nodeList = children_.get(name);
		if (nodeList == null)
		{
			// Not necessarily an exception if there are no child nodess
			System.err.println("getNode: I am " + name_
					+ " and I have no child node by the name " + name);
			System.err.print("getNode: My children are:");
			for (String childName : children_.keySet())
				System.err.print(" " + childName);
			System.err.println();
			return null;
		}
		return nodeList.get(index);
	}
	
	/**
	 * Returns the first node with the specified name. Equivalent to
	 * getNode(name, 0).
	 * 
	 * @param name
	 * @return
	 */
	public MaeNode getNode(String name)
	{
		return getNode(name, 0);
	}

	/**
	 * Is this node a data array? If it is indeed a data array, it doesn't have
	 * any child nodes.
	 * 
	 * @return True if this node is a data array
	 */
	public boolean isDataArray()
	{
		return isDataArray_;
	}

	/**
	 * This MaeNode can have zero or more child MaeNodes. Returns a list of
	 * child node names.
	 * 
	 * @return The names of the child nodes, if any
	 */
	public Set<String> getChildNodeNames()
	{
		return children_.keySet();
	}

	/**
	 * Returns the name of this node. Special names include:
	 * anonymous: The first block in a .mae file that includes the version
	 *   of the format
	 * HEAD: Node that encloses all the other nodes and is not actually part
	 *   of the mae file
	 * 
	 * @return The name of this node
	 */
	public String getName()
	{
		return name_;
	}

	/**
	 * Returns a list of data headings. It's a list because you might care about
	 * their ordering (though this is only relevant for cosmetic purposes).s
	 * 
	 * @return The headings (field names) for the data record(s)
	 */
	public List<String> getHeadings()
	{
		return headings_;
	}

	/**
	 * Returns the element in the specified row (index) in the specified column.
	 * 
	 * e.g.: getNode("f_m_ct", 1).getString("s_ffio_ct_type", 1); firstAtomX =
	 * getNode("f_m_ct", 1).getNode("m_atom", 1).getInt("r_m_x_coord", 1);
	 * 
	 * @param columnName
	 * @param index
	 * @return The element you asked for, or null if it couldn't find it
	 */
	public String getString(String columnName, int index)
	{
		ArrayList<String> column = dataColumns_.get(columnName);
		if (column == null)
			return null;
		else
			return column.get(index);
	}

	/**
	 * Makes sure that the backing data structures for the specified column
	 * exist, so subsequent calls to setString and whatnot will work.
	 * 
	 * @param columnName
	 */
	private void ensureColumnExists(String columnName)
	{
		ArrayList<String> column = dataColumns_.get(columnName);
		if (column == null)
		{
			// If the column doesn't exist, add it
			column = new ArrayList<String>();
			dataColumns_.put(columnName, column);
			// XXX: This is kind of dumb. We have two lists of headings, one in dataColumns_
			// and one in headings_, that we have to keep in sync.
			if(headings_.contains(columnName) == false)
				headings_.add(columnName);
			// If a column didn't exist, that means we are in a constructor
			// and there are no complete data records yet
		}
	}

	/**
	 * Sets a data element value. Requires the column in question to either not
	 * exist or already be big enough, because I am too lazy to write the code
	 * to resize all the other columns. In the case where the column didn't
	 * already exist, we're probably being called from a constructor and so it's
	 * OK to make the column here.
	 * 
	 * @param columnName
	 *            Heading (field name) of the column
	 * @param index
	 *            Row number in the column
	 * @param value
	 *            Value to set it to
	 */
	public void setString(String columnName, int index, String value)
	{
		ensureColumnExists(columnName);
		ArrayList<String> column = dataColumns_.get(columnName);
		if (column.size() == 0)
			column.add(value);
		else if (column.size() > index)
			column.set(index, value);
		else
		{
			System.err.println("Node.setString: Index out of bounds! How's that for a cryptic error?");
			throw new IndexOutOfBoundsException();
		}

	}

	/**
	 * Sets the value of the first row in a column. This is mainly useful if
	 * this MaeNode is *not* a data array and thus there is only one data
	 * record.
	 * 
	 * @param columnName
	 * @param value
	 */
	public void setString(String columnName, String value)
	{
		setString(columnName, 0, value);
	}

	/**
	 * Like getString, but does the parseDouble for you.
	 * 
	 * @param columnName
	 * @param index
	 * @see getString
	 * @return The value you asked for
	 */
	public double getDouble(String columnName, int index)
	{
		return Double.parseDouble(getString(columnName, index));
	}

	/**
	 * Like getString, but does the parseInt for you.
	 * 
	 * @param columnName
	 * @param index
	 * @see getString
	 * @return The value you asked for
	 */
	public int getInteger(String columnName, int index)
	{
		return Integer.parseInt(getString(columnName, index));
	}

	/**
	 * Links a node to this one as a child.
	 * 
	 * @param node
	 */
	public void addChildNode(MaeNode node)
	{
		if (isDataArray_)
		{
			System.err
					.println("addChildNode: This node is a data array, which means a child node cannot be added to it.");
			return;
		}

		String nodeName = node.getName();
		if (nodeName == null) // Nodes with no name are not really nodes
			return;
		// Only the head node can have two child nodes of the same name, so
		// we disallow it here
		if (children_.containsKey(nodeName) == false)
			children_.put(nodeName, new ArrayList<MaeNode>());
		ArrayList<MaeNode> nodes = children_.get(nodeName);
		nodes.add(node);
		// System.err.println("addChildNode: Added node " + nodeName + " to " +
		// name_);
	}

	/**
	 * Empty constructor.
	 */
	private MaeNode()
	{
	}

	/**
	 * General constructor for MaeNode.
	 * 
	 * @param name
	 *            Name of the node (e.g. "f_m_ct")
	 * @param isDataArray
	 *            Specifies whether this node is a data array. If so, it cannot
	 *            have any child nodes.
	 * @param headings
	 *            List of headings for the data in this node
	 */
	public MaeNode(String name, boolean isDataArray, String... headings)
	{
		name_ = name;
		for (String heading : headings)
			headings_.add(heading);
		isDataArray_ = isDataArray;
		if (isDataArray == false)
		{
			// Got to have data, even if just placeholders
			for (String heading : headings_)
				setString(heading, heading + "_PLACEHOLDER");
			numDataRecords_ = 1;
		}
	}

	/**
	 * Reads a Mae format (.mae/.cms) file and recursively reconstructs the tree of nodes
	 * in that file.
	 * 
	 * @param filename
	 * @return
	 * @throws IOException
	 */
	public static MaeNode fromCMSFile(String filename) throws IOException
	{
		StreamTokenizer input = new StreamTokenizer(new BufferedReader(
				new FileReader(filename)));

		// Set up StreamTokenizer, one of the most gimpy classes in Java
		input.eolIsSignificant(false);
		input.wordChars(0x21, 255);
		input.wordChars('_', '_');
		input.ordinaryChars('%', '%');
		input.wordChars('%', '%');
		input.ordinaryChars('\'', '\'');
		input.wordChars('\'', '\'');
		input.ordinaryChars('.', '.');
		input.ordinaryChars(':', ':');
		input.ordinaryChars('-', '-');
		input.wordChars('-', '-');
		input.ordinaryChars('0', '9');
		input.quoteChar('"');
		input.wordChars('0', '9');
		input.wordChars('.', '.');
		input.wordChars(':', ':');
		input.commentChar('#');

		MaeNode head = new MaeNode();
		head.name_ = "HEAD"; // Magic name for root of the tree of nodess

		// Keep reading data form the stream until EOF
		while (input.ttype != StreamTokenizer.TT_EOF)
		{
			MaeNode node = new MaeNode(input);
			head.addChildNode(node);
			// System.err.println("=========== (toplevel)");
		}

		return head;
	}

	
	
	/**
	 * Searches for the f_m_ct child node with the specified field matching specified value.
	 * This is useful for finding solute, solvent, etc blocks.
	 * 
	 * @returns The first f_m_ct node matching the criterion, or null if none were found.
	 */
	private MaeNode getCTNodeByField(String field, String value) throws MaeNodeException
	{
		if (!getName().equals("HEAD"))
			throw new MaeNodeException(
					"You must specify a HEAD node, because no other node type contains f_m_ct nodes");

		ArrayList<MaeNode> nodes = children_.get("f_m_ct");
		for (MaeNode node : nodes)
		{
			String theValue = node.getString(field, 0);
			if (theValue.equals(value))
				return node;
		}

		return null;
	}
	
	public MaeNode getCTNodeByType(String type) throws MaeNodeException
	{
		return getCTNodeByField("s_ffio_ct_type", type);
	}
	
	public MaeNode getCTNodeByTitle(String title) throws MaeNodeException
	{
		return getCTNodeByField("s_m_title", title);
	}
	
	/**
	 * Maestro is completely worthless at exporting its own format to PDB when
	 * the structure is not valid by its standards (e.g. too many bonds).
	 * 
	 * Hence this function.
	 * 
	 * @param out Stream to output to, like a file or System.out
	 * @throws MaeNodeException If the MaeNode is screwed up
	 */
	public void soluteOnlyToStreamAsPDB(PrintStream out) throws MaeNodeException
	{
		MaeNode solute = getCTNodeByType("solute");
		MaeNode solute_m_atom = solute.getNode("m_atom");
		
		if (solute_m_atom == null)
			throw new MaeNodeException("soluteOnlyTotreamAsPDB: m_atom is missing! Nuts!");
		
		out.println("REMARK 888 From MaeTool by Tom Joseph <thomas.joseph@mssm.edu>");
		out.println("REMARK 888 Cool, eh?");
		
		int numAtoms = solute_m_atom.numDataRecords();
				
		for(int atom = 0; atom < numAtoms; atom++)
		{
			String name = solute_m_atom.getString("s_m_pdb_atom_name", atom).trim();
			String resName = solute_m_atom.getString("s_m_pdb_residue_name", atom).trim();
			//char chain = solute_m_atom.getString("s_m_chain_name", atom).trim().charAt(0);
			char atomNamePrefix = ' ';
			if(name.length() > 3)
			{
				atomNamePrefix = name.charAt(0);
				name = name.substring(1);
			}
			
			char chain = 'X';
			int resID = solute_m_atom.getInteger("i_m_residue_number", atom);
			double x = solute_m_atom.getDouble("r_m_x_coord", atom);
			double y = solute_m_atom.getDouble("r_m_y_coord", atom);
			double z = solute_m_atom.getDouble("r_m_z_coord", atom);
			double occupancy = 0.0;
			double tempFactor = 0.0;
			out.printf(
			//          11111111112222222222333333333344444444445555555555666666666677777777778
			// 123456789012345678901234567890123456789012345678901234d6789012345678901234567890
			//"ATOM  %5d--_%4s-_%4s%c%4d_-   %8.3f---%8.3f---%8.3f---%6.2--%6.2                ",
			  "ATOM  %5d %c%-3s %-3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f              \n",
			  atom + 1,atomNamePrefix, name, resName, chain, resID, x, y, z, occupancy, tempFactor);
		}
		out.println("END");
	}

	/**
	 * Writes the contents of this node as a PSF in the X-PLOR format. (NAMD
	 * requires the X-PLOR format.) For this to work, this node must contain the
	 * following blocks: m_atom m_bond
	 * 
	 * @param out Stream to write to
	 * @throws MaeNodeException If the MaeNode is screwed up
	 */
	public void soluteOnlyToStreamAsPSF(PrintStream out) throws MaeNodeException
	{
		MaeNode solute = getCTNodeByType("solute");
		MaeNode solute_m_atom = solute.getNode("m_atom");
		MaeNode solute_m_bond = solute.getNode("m_bond");

		if (solute_m_atom == null || solute_m_bond == null)
			throw new MaeNodeException("soluteOnlyToStreamAsPSF: m_atom and/or m_bond are missing!");

		// We'll just examine the PDB name column
		ArrayList<String> atomNames = solute_m_atom.getColumn("s_m_pdb_atom_name");
		if (atomNames == null)
			throw new MaeNodeException(
					"writeMinimalPSF: This m_atom block doesn't seem to have the s_m_pdb_atom_name field.  I need that.");

		// PSF global header
		out.println("PSF CMAP\n");
		out.println("       2 !NTITLE");
		out.println(" REMARKS Generated by MaeTool, written by Tom Joseph");
		out.println(" REMARKS <thomas.joseph@mssm.edu>\n");

		// Want to do the following sections:
		// NBOND
		// NTHETA: angles
		// NPHI: proper dihedrals
		// NIMPHI: improper dihedrals
		// NCRTERM: cross-terms

		int numAtoms = solute_m_atom.numDataRecords();
		out.printf("%1$8d !NATOM\n", numAtoms);

		for (int atom = 0; atom < numAtoms; atom++)
		{
			String segID = solute_m_atom.getString("s_m_chain_name", atom).trim();
			String resID = solute_m_atom.getString("i_m_residue_number", atom).trim();
			String resName = solute_m_atom.getString("s_m_pdb_residue_name", atom).trim();
			String atomName = solute_m_atom.getString("s_m_pdb_atom_name", atom).trim();

			// atomid segname resid resname atomname atomtype charge mass 0.0
			out.printf(
					"%1$8d %2$-4s %3$-4d %4$-4s %5$-4s %6$-4s %7$12.6f %8$10.4f\n",
					atom + 1, segID, Integer.parseInt(resID), resName,
					atomName, atomName, 0.0, 1.0);
		}

		// TODO: solvent and ions - must offset atom IDs
		// Order of atoms within full_system block is solute, solvent, ions

		// NBOND
		out.println();
		int numBonds = solute_m_bond.numDataRecords();
		out.printf("%1$8d !NBOND\n", numBonds);

		for (int bond = 0; bond < numBonds; bond++)
		{
			int from = solute_m_bond.getInteger("i_m_from", bond);
			int to = solute_m_bond.getInteger("i_m_to", bond);
			out.printf("%1$8d%2$8d", from, to);
			// Four pairs of atoms per line
			if (bond > 0 && (bond % 4 == 3))
				out.println();
		}

		// NTHETA (angles)
		out.println();
		// This has to come from the ffio_ff block
		MaeNode solute_ffio_angles = solute.getNode("ffio_ff").getNode(
				"ffio_angles");
		int numAngles = solute_ffio_angles.numDataRecords();
		out.printf("%1$8d !NTHETA\n", numAngles);
		// Print the angles...
		for (int angle = 0; angle < numAngles; angle++)
		{
			int ai = solute_ffio_angles.getInteger("i_ffio_ai", angle);
			int aj = solute_ffio_angles.getInteger("i_ffio_aj", angle);
			int ak = solute_ffio_angles.getInteger("i_ffio_ak", angle);
			out.printf("%1$8d%2$8d%3$8d", ai, aj, ak);
			// Three triples of atotms per line
			if (angle > 0 && (angle % 3 == 2))
				out.println();
		}
		
		// NPHI and NIMPHI
		// Mae files have both proper and improper dihedrals in the same block.
		// So we need to scan the block and extract each kind.
		ArrayList<int[]> propers = new ArrayList<int[]>(), 
			impropers = new ArrayList<int[]>();
		MaeNode solute_ffio_dihedrals = solute.getNode("ffio_ff").getNode("ffio_dihedrals");
		int numDihedrals = solute_ffio_dihedrals.numDataRecords();
		
		for (int dihedral = 0; dihedral < numDihedrals; dihedral++)
		{
			String type = solute_ffio_dihedrals.getString("s_ffio_funct", dihedral);
			int ai = solute_ffio_dihedrals.getInteger("i_ffio_ai", dihedral);
			int aj = solute_ffio_dihedrals.getInteger("i_ffio_aj", dihedral);
			int ak = solute_ffio_dihedrals.getInteger("i_ffio_ak", dihedral);
			int al = solute_ffio_dihedrals.getInteger("i_ffio_al", dihedral);
			int fourAtoms[] = {ai, aj, ak, al};
			if(type.equals("Improper_Trig"))
				impropers.add(fourAtoms);
			else if(type.equals("Proper_Trig"))
				propers.add(fourAtoms);				
		}
		
		out.println();
		int numPropers = propers.size();
		out.printf("%1$8d !NPHI\n", numPropers);
		for (int proper = 0; proper < numPropers; proper++)
		{
			int[] record = propers.get(proper);
			out.printf("%1$8d%2$8d%3$8d%4$8d", record[0], record[1], record[2], record[3]);
			// Two quadruples per line
			if (proper > 0 && (proper % 2 == 1))
				out.println();
		}

		out.println();
		int numImpropers = impropers.size();
		out.printf("%1$8d !NIMPHI\n", numPropers);
		for (int improper = 0; improper < numImpropers; improper++)
		{
			int[] record = impropers.get(improper);
			out.printf("%1$8d%2$8d%3$8d%4$8d", record[0], record[1], record[2], record[3]);
			// Two quadruples per line
			if (improper > 0 && (improper % 2 == 1))
				out.println();
		}
		
		// No cross-terms; I don't even know what those are
		out.printf("\n%1$8d !NCRTERM\n\n", 0);

	}

	/**
	 * Prints this Node and its children to a stream. This is meant to generate
	 * a file that Desmond/Maestro can read.
	 * 
	 * @param out
	 *            A stream where the output should go, such as System.out
	 */
	public void toStream(PrintStream out)
	{
		// The HEAD node doesn't exist in the file; it's only there to contain
		// all the other nodes
		if (name_.equals("HEAD") == false)
		{
			// anonymous node is the first node in the file and
			// contains the format version
			if (name_.equals("anonymous") == false)
			{
				out.print(indent_ + name_);
				if (isDataArray_)
					out.print("[" + numDataRecords_ + "]");
				out.print(" ");
			}
			out.println("{");
			increaseIndent();

			// Print headings
			for (String heading : headings_)
				out.println(indent_ + heading);
			out.print(indent_ + ":::");
			if (isDataArray_)
				out.println();
		}

		// There is always at least one data record, even if the node is not
		// a data array
		for (int i = 0; i < numDataRecords_; i++)
		{
			out.print(indent_);
			if (isDataArray_)
				out.print(i + 1); // Array index

			for (String heading : headings_)
			{
				if (isDataArray_)
					out.print(" ");
				else
					out.print("\n" + indent_);
				String s = getString(heading, i);
				if (s.length() == 0 || s.indexOf(' ') >= 0)
					out.print("\"" + s + "\"");
				else
					out.print(s);
			}
			out.println();
		}

		// Data arrays have an extra ::: at the end of the data
		if (isDataArray_)
			out.println(indent_ + ":::");

		// Iterate over each name and output all nodes of that name (there can
		// be more
		// than one node with the same name, e.g. f_m_ct)
		for (String nodeName : children_.keySet())
		{
			ArrayList<MaeNode> nodes = children_.get(nodeName);
			for (MaeNode node : nodes)
				node.toStream(out);
		}

		// Don't print the HEAD node - it is just an internal convenience
		if (name_.equals("HEAD") == false)
		{
			decreaseIndent();
			out.println(indent_ + "}");
		}
	}

	/**
	 * Creates this node from an input .mae-formatted stream. Presumably the 
	 * StreamTokenizer was created from a file. We do things this way to abstract 
	 * away the file manipulation stuff from the constructor, so that this
	 * constructor can easily recurse.
	 * 
	 * @param input
	 * @throws IOException
	 */
	public MaeNode(StreamTokenizer input) throws IOException
	{
		// Read in the data table and recurse on child nodes

		// To begin, we are starting at the beginning of a node.
		// The first node in a file has no name.
		// Grab a token. If it is a '{', the node has no name.
		input.nextToken();

		//System.err.println(input.lineno() + ": The first token in this node is: " + input.sval);
		// Blank line(s) at end of file will cause input.sval to be null, I
		// guess
		if (input.sval == null || input.ttype == StreamTokenizer.TT_EOF
				|| input.sval.equals("}"))
		{
			// There is not actually a node here
			name_ = null;
			// System.err.println("EXITING RECURSION at beginning of constructor: "
			// + input.sval);
			return;
		}

		if (input.sval.equals("{") == false)
		{
			name_ = input.sval;
			// If the last token was not a '{', the next token is, so eat it.
			input.nextToken();
			// System.err.println("Should have eaten {: " + input.sval);
			assert (input.sval.equals("{"));
		}

		// At this point we are positioned after the opening brace

		// TODO:
		// If the node has no name, check if the version is what we are
		// expecting
		// (that is, s_m_m2io_version is 2.0.0)

		// If the name of this node has a [] suffix, it's an array and nothing
		// else.
		numDataRecords_ = 1;
		if (name_.endsWith("]"))
		{
			int a = name_.indexOf('[') + 1;
			int b = name_.indexOf(']');
			numDataRecords_ = Integer.parseInt(name_.substring(a, b));
			// System.err.println("Reading " + numDataRecords_ +
			// " data records.");
			isDataArray_ = true;
			name_ = name_.substring(0, a - 1);
		}
		// System.err.println("This node is named " + name_);
		// Otherwise, it's a node with a single keys/values section and zero or
		// more child nodes.
		// Read column headings
		int numHeadings = 0;
		input.nextToken(); // Get the first heading value
		while (input.sval.equals(":::") == false)
		{
			// System.err.println("Heading: " + input.sval);
			dataColumns_.put(input.sval, new ArrayList<String>());
			headings_.add(input.sval);
			numHeadings++;
			input.nextToken();
		}

		// Read however many lines of values we are supposed to. If this is an
		// array, there is
		// a ':::' token at the end that we need to eat
		//System.err.println("Reading " + numDataRecords_ + " records of "
		//		+ numHeadings + " values each.");
		for (int valueSet = 0; valueSet < numDataRecords_; valueSet++)
		{
			// Eat array index number
			if (isDataArray_)
				input.nextToken();

			for (int heading = 0; heading < numHeadings; heading++)
			{
				input.nextToken();
				String headingName = headings_.get(heading);
				dataColumns_.get(headingName).add(input.sval);
			}

			// System.err.println("-");
		}

		if (isDataArray_)
			input.nextToken(); // Eat ":::"

		input.nextToken(); // Get either } or next node name

		// Now there will be zero or more child nodes
		while (input.sval.equals("}") == false)
		{
			// System.err.println("Starting new node on token: " + input.sval);
			input.pushBack();
			MaeNode node = new MaeNode(input);
			addChildNode(node);
			input.nextToken(); // Eat closing brace
			// System.err.println("-------------");
		}

		// System.err.println("EXITING RECURSION at end of constructor: " +
		// input.sval);
	}

	/**
	 * Returns the number of data records.
	 * 
	 * @return
	 */
	public int numDataRecords()
	{
		return numDataRecords_;
	}

	/**
	 * Adds a data record to this node. You need to know the names of the
	 * headings to use this. You can use getHeadings for that.
	 * 
	 * @param record
	 *            A map from keys to values for data you want to add.
	 * @see getHeadings
	 */
	public void addDataRecord(Map<String, String> record)
	{
		for (String heading : headings_)
		{
			if (record.containsKey(heading) == false)
			{
				System.err
						.println("addDataRecord: The record I was passed doesn't have key "
								+ heading + ", but it really should.");
				System.exit(1);
			}

			// Add the element to this column
			ensureColumnExists(heading);
			dataColumns_.get(heading).add(record.get(heading));
		}
		numDataRecords_++;
	}
	
	/**
	 * Deletes a data record, given its index; counting starts at zero.
	 * 
	 * @param index
	 * @throws MaeNodeException If there isn't a record with that index.
	 */
	public void deleteDataRecord(int index) throws IndexOutOfBoundsException
	{
		// Does this record exist?
		if(index >= numDataRecords())
		{
			throw new IndexOutOfBoundsException("deleteDataRecord: Tried to delete record " + index  
					+ " but there are only " + numDataRecords() + " records in this node");
		}
		
		for(String heading : headings_)
		{
			dataColumns_.get(heading).remove(index);
		}
		
		numDataRecords_--;
	}

	/**
	 * Returns a reference to an entire column. It's more efficient to use this
	 * instead of the other get*() if you are just interested in the contents of a single
	 * column.
	 * 
	 * @param string
	 * @return Reference to the column
	 */
	public ArrayList<String> getColumn(String string)
	{
		return dataColumns_.get(string);
	}
	
	/**
	 * Sets up an FEP calculation using Desmond.
	 * 
	 * This node must have a single solute CT block. You must specify
	 * 
	 * @param perturbed A CMS or whatever containing a solute CT block that
	 *     will serve as the perturbed version.
	 * @param phase_out Atom IDs in the original solute that will become 
	 *     dummy atoms during the perturbation
	 * @param phase_in Atom IDs in the perturbed solute that will start as 
	 *     dummy and become real during the perturbation
	 */
	public void setupForAlchemicalFEP(MaeNode perturbedCMS)
		throws MaeNodeException
	{
		// assert that we are the HEAD node
		if(getName().equals("HEAD") == false
				|| perturbedCMS.getName().equals("HEAD") == false)
			throw new MaeNodeException("This method only works on two entire CMS files.");
		
		// The ffio_ff blocks should be within their respective ct blocks
		MaeNode original = getCTNodeByType("solute");
		MaeNode perturbed = perturbedCMS.getCTNodeByType("solute");
		if(original.getNode("ffio_ff") == null || perturbed.getNode("ffio_ff") == null)
			throw new MaeNodeException("Missing one or both ffio_ff blocks.");
				
		// Add s_fepio_name ct-level fields
		// Apparently this can be a friendly, descriptive name
		original.setString("s_fepio_name", "Original");
		perturbed.setString("s_fepio_name", "Perturbed");
		
		// Add i_fepio_stage ct-level fields
		original.setString("i_fepio_stage", "1");
		perturbed.setString("i_fepio_stage", "2");
		
		// The perturbed ct block must have an fepio_fep block that contains
		// the mappings from atoms, bonds, angles, etc from original to
		// perturbed structures.
		MaeNode fepio_fep = new MaeNode("fepio_fep", false, "s_fepio_name", "i_fepio_stage");
		// Point to the *original* ct block
		fepio_fep.setString("s_fepio_name", "Original");
		fepio_fep.setString("i_fepio_stage", "1");
		
		// Now we need to construct the mappings from original to perturbed.
		// Our strategy will assume that the perturbation won't cause a shift in residue
		// numbers (i_m_residue_number field in m_atom[]).
		// All atoms/bonds/whatever should be the same between original and perturbed, except
		// in the perturbed region, which we will just mutate to/from dummy atoms.
		// We will compare s_m_pdb_residue_name, s_m_pdb_atom_name, and r_m_[xyz]_coord fields
		// between each atom in each residue.
		
		// Make list of residues that are different between original and perturbed.
		// Residues with the same i_m_residue_number are supposed to be the same.
		MaeNode origAtoms = original.getNode("m_atom");
		MaeNode pertAtoms = perturbed.getNode("m_atom");
		ArrayList<Integer> disappearingAtoms = new ArrayList<Integer>();
		ArrayList<Integer> appearingAtoms = new ArrayList<Integer>();
		
		int index = 0, origResidueStart = 0, origResidueEnd = -1;
		int pertResidueStart = 0, pertResidueEnd = -1;
		// Compare the residues. If they aren't the same, mark them as "transition" residues
		HashMap<Integer, Integer> singleResidueAtomMap, atomMap = new HashMap<Integer, Integer>();

		// Iterate over atoms
		// XXX: Assumes same number of residues and 1-to-1 mapping between residues with the
		// same residue ID!!! This is pretty lame but it makes things easier, and covers my use case.
		while(origResidueStart >= 0 && pertResidueStart >= 0)
		{
			// Establish bounds of current residue in both original and perturbed m_atom blocks
			// The first residue must start at the first atom
			origResidueStart = origResidueEnd + 1;
			origResidueEnd = findResidueEnd(origAtoms, origResidueStart);
			pertResidueStart = pertResidueEnd + 1;
			pertResidueEnd = findResidueEnd(pertAtoms, pertResidueStart);

			// Map the atoms within each residue to each other
			singleResidueAtomMap = mapAtomsBetweenResidues(origAtoms, origResidueStart, origResidueEnd,
					pertAtoms, pertResidueStart, pertResidueEnd);

			// If the construction of the map failed, the two residues must be different,
			// so they must be part of the alchemical mutation
			if(singleResidueAtomMap == null)
			{
				// Because of the way the atom map numbering works, we must save the atom IDs
				// and construct this part of the atom map later.
				for(int ai = origResidueStart; ai <= origResidueEnd; ai++)
					disappearingAtoms.add(ai);

				for(int aj = pertResidueStart; aj <= pertResidueEnd; aj++)
					appearingAtoms.add(aj);
			} else
			{
				// Save this atom mapping so we can add it to fepio_atommaps later
				atomMap.putAll(singleResidueAtomMap);
			}
		} // End iteration over atoms
		
		// Now we do the atom maps for the perturbed atoms.
		// Each atom disappearing in the original structure gets mapped to -1.
		for(Integer i : disappearingAtoms)
			atomMap.put(i, -1);

		// Each atom appearing in the perturbed structure gets mapped to -(#atoms_in_original + i)
		int numAtomsInOriginal = origAtoms.numDataRecords();
		for(Integer i : appearingAtoms)
			atomMap.put(-(numAtomsInOriginal + i), i);
		
		// Dump atomMap to fepio_atommaps
		MaeNode fepio_atommaps = new MaeNode("fepio_atommaps", true, "i_fepio_ai", "i_fepio_aj");
		HashMap<String, String> tmpMap = new HashMap<String, String>();
		for (Integer ai : atomMap.keySet())
		{
			tmpMap.put("i_fepio_ai", ai.toString());
			tmpMap.put("i_fepio_aj", atomMap.get(ai).toString());
			fepio_atommaps.addDataRecord(tmpMap);
		}
		
		// Now we need to do fepio_bondmaps.
		//
		// i_fepio_ti: indexes bond in original CT ffio_bonds block
		// i_fepio_tj: the corresponding bond in perturbed CT ffio_bonds block
		// i_fepio_ai: First atom in this bond in the original CT. Can be negative as in atommaps
		// i_fepio_aj: Second atom in this bond in the original CT.
		//
		// While we assumed original and perturbed have the same number of *residues*, they may
		// not have the same number of *atoms* and therefore ffio_bonds/m_bond for the solute
		// may be different. There will probably be a shift at some point to accommodate an
		// insertion of atoms in one of them.
		//
		// So, for each bond potential in original, we must determine which two atoms in perturbed it
		// refers to (accomplished using atom map we constructed above) and using this information,
		// find the bond potential in perturbed that the original bond potential must map to.
		
		MaeNode origBonds = original.getNode("ffio_ff").getNode("ffio_bonds");
		MaeNode pertBonds = perturbed.getNode("ffio_ff").getNode("ffio_bonds");
		HashMap<Integer, Integer> bondMap = new HashMap<Integer, Integer>();
				
		int numOrigBonds = origBonds.numDataRecords(), startFrom = 0;
		int numPertBonds = pertBonds.numDataRecords();
		for(int origBondNum = 0; origBondNum < numOrigBonds; origBondNum++)
		{
			// Get original atoms and the perturbed atoms they correspond tos
			int origAi = origBonds.getInteger("i_ffio_ai", origBondNum);
			int origAj = origBonds.getInteger("i_ffio_aj", origBondNum);
			int pertAi = atomMap.get(origAi);
			int pertAj = atomMap.get(origAj);
			
			// If any of the original atoms are dummy atoms, we need to say
			// that this bond does not exist in the perturbed state by specifying
			// a negative bond potential number (I will use -1 because the docs
			// do not provide specific guidance otherwise).
			if(origAi < 0 || origAj < 0)
			{
				bondMap.put(origBondNum, -1);
				continue;
			}
			
			// Find the bond from pertAi to pertAj in pertBonds
			// We continue to advance where we start looking as this is just a linear
			// search and we expect this to reduce the average amount of time spent
			// looking for the correct bond.
			boolean foundPertBond = false;
			for(int pertBondNum = startFrom + 1; 
				pertBondNum != startFrom; 
				pertBondNum = (pertBondNum + 1) % numPertBonds)
			{
				if(pertBonds.getInteger("i_ffio_ai", pertBondNum) == pertAi
						&& pertBonds.getInteger("i_ffio_aj", pertBondNum) == pertAj)
				{
					bondMap.put(origBondNum, pertBondNum);
					startFrom = pertBondNum;
					foundPertBond = true;
					// Technically we don't need this break as the above assignment will
					// trigger the termination condition of the for loop
					break;
				}
			}
			
			// If we didn't find the corresponding bond, something is corrupt
			if(foundPertBond == false)
				throw new MaeNodeException("Couldn't find corresponding perturbed bond for original bond " + origBondNum);
			
		} // End iteration over original bonds
		
		// Now, fill out fepio_bondmaps
		MaeNode fepio_bondmaps = new MaeNode("fepio_bondmaps", true,
				"i_fepio_ti", "i_fepio_tj", "i_fepio_ai", "i_fepio_aj");
		// TODO: actually fill out bond maps

	}
	
	/**
	 * Returns a mapping if the residue in m_atom block atoms1, starting at index1, is the "same"
	 * as the residue starting at index2 in atoms2. This isn't particularly robust but it
	 * should be good enough for FEP stuff...
	 * 
	 * The mapping can be used for fepio_atommaps.
	 * 
	 * XXX: Assumes the atoms are ordered the same in both residues!
	 * 
	 * @param atoms1 m_atom block
	 * @param start1 Index within atoms1 at which the residue starts
	 * @param end1
	 * @param atoms2
	 * @param start2
	 * @param end2
	 * @return
	 */
	private static HashMap<Integer, Integer> mapAtomsBetweenResidues(MaeNode atoms1, int start1, int end1, 
			MaeNode atoms2, int start2, int end2)
	{
		int length1 = end1 - start1 + 1, length2 = end2 - start2 + 1;
		// If the residues have different numbers of atoms they must be different
		if(length1 != length2)
			return null;
		
		HashMap<Integer, Integer> atomMap = new HashMap<Integer, Integer>();
		
		// Iterate through each atom of the first residue and see if there is a matching
		// atom in the other residue.
		// (Yes, this is O(n^2) but n should be small.)
		for(int i = start1; i <= end1; i++)
		{
			boolean foundIt = false;
			String resName = atoms1.getString("s_m_pdb_residue_name", i).trim();
			String atomName = atoms1.getString("s_m_pdb_atom_name", i).trim();
			double x = atoms1.getDouble("r_m_x_coord", i);
			double y = atoms1.getDouble("r_m_y_coord", i);
			double z = atoms1.getDouble("r_m_z_coord", i);
			// We always start on the "equivalent" atom index and wrap around so the common
			// case of identical atom order within identical residues is optimized
			int offset = i - start1; 
			for(int j = 0; j < length2; j++)
			{
				int idx = (j + offset) % length2 + start2;
				String resName2 = atoms2.getString("s_m_pdb_residue_name", idx).trim();
				String atomName2 = atoms2.getString("s_m_pdb_atom_name", idx).trim();
				double x2 = atoms2.getDouble("r_m_x_coord", idx);
				double y2 = atoms2.getDouble("r_m_y_coord", idx);
				double z2 = atoms2.getDouble("r_m_z_coord", idx);
				// Compare s_m_pdb_residue_name, s_m_pdb_atom_name, and r_m_[xyz]_coord fields
				if( resName.equals(resName2) &&	atomName.equals(atomName2)
					&& x == x2 && y == y2 && z == z2 )
				{
					foundIt = true;
					System.err.println("mapAtomsBetweenResidues: " + i + " <-> " + idx);
					atomMap.put(i, idx);
					break;
				} else
				{
					System.err.println("mapAtomsBetweenResidues: " + i + " does not map to " + idx);
				}
			}
			
			// If we couldn't find a matching atom, give up
			if(foundIt == false)
				return null;			
		}
				
		return atomMap;
	}
	
	// XXX Maybe don't need this?
	private static int findResidueStart(MaeNode m_atom, int startingFrom) throws MaeNodeException
	{
		// Start at startingFrom and move backwards until the i_m_residue_number changes
		int currentResidueNumber = m_atom.getInteger("i_m_residue_number", startingFrom);
		int index = startingFrom, lastResidueNumber = currentResidueNumber;
		while(lastResidueNumber == currentResidueNumber && index >= 0)
		{
			index--;
			lastResidueNumber = m_atom.getInteger("i_m_residue_number", index);
		}
		return index + 1;
	}

	/**
	 * Finds the next ending of a residue with respect to a specified index.
	 * 
	 * @param m_atom Atom block to look at
	 * @param startingFrom Where in the atom block to start looking
	 * @return Index of the last atom in this residue
	 * @throws MaeNodeException
	 */
	private static int findResidueEnd(MaeNode m_atom, int startingFrom) throws MaeNodeException
	{
		int numAtoms = m_atom.numDataRecords();
		if(startingFrom > numAtoms)
			return -1;
		// Start at startingFrom and move forward until the i_m_residue_number changes
		int currentResidueNumber = m_atom.getInteger("i_m_residue_number", startingFrom);
		int index = startingFrom, lastResidueNumber = currentResidueNumber;
		while(lastResidueNumber == currentResidueNumber && index < (numAtoms - 1))
		{
			index++;
			lastResidueNumber = m_atom.getInteger("i_m_residue_number", index);
		}
		return index - 1;
	}

	/**
	 * Sets all the fields of a particular data record at once.
	 * 
	 * @param record A map of key-value pairs. If the not all headings are present here
	 *   keys, the value will default to "<>", which means blank
	 * @param index Index of the data record to set
	 */
	public void setDataRecord(HashMap<String, String> record, int index)
	{
		for (String heading : headings_)
		{
			if (record.containsKey(heading) == false)
			{
				System.err.println("setDataRecord: The record I was passed doesn't have key "
								+ heading + ", but it really should.");
				System.err.println("setDataRecord: Defaulting to <>");
				record.put(heading, "<>");
			}

			// Add the element to this column
			ensureColumnExists(heading);
			//dataColumns_.get(heading).add(record.get(heading));
			setString(heading, index, record.get(heading));
		}
	}

	/**
	 * Given an atom ID (counted starting at 1), delete that atom
	 * and all bonds associated with it from m_bond and m_atom. Does not touch
	 * ffio_ff block, so this method effectively breaks it.
	 * 
	 * @param atomID The atom ID to delete (counted starting at 1)
	 * @throws MaeNodeException If this node doesn't have both m_bond and m_atom child nodes.
	 * 	In general this node should be a CT block
	 */
	public void deleteAtom(int atomID) throws MaeNodeException
	{
		MaeNode m_bond = getNode("m_bond");
		MaeNode m_atom = getNode("m_atom");
		
		if(m_bond == null || m_atom == null)
		{
			throw new MaeNodeException("deleteAtom: The specified node doesn't have m_atom and/or m_bond blocks.");
		}
		
		// Find and nuke any bonds involving this atom
		int numBonds = m_bond.numDataRecords();
		for(int bond = 0; bond < numBonds; bond++)
		{
			int fromIndex = m_bond.getInteger("i_m_from", bond);
			int toIndex = m_bond.getInteger("i_m_to", bond);

			if(fromIndex == atomID || toIndex == atomID)
			{
				m_bond.deleteDataRecord(bond);
				// The rest of the records will shift up since they are stored
				// in ArrayLists
				numBonds--;
				
				// Since we shifted the rest of the bond records up, there is now
				// a different bond with this bond index, so we adjust our iterating
				// index variable so we don't skip it
				bond--;
			} else
			{
				// Shift atom indices referenced by bonds as necessary
				// If the atom index referenced is greater than the atom index
				// we are deleting, subtract one from it
				if(fromIndex > atomID)
					m_bond.setInteger("i_m_from", bond, fromIndex - 1);
				
				if(toIndex > atomID)
					m_bond.setInteger("i_m_to", bond, toIndex - 1);
			}
		}
		
		// Delete the atom itself
		m_atom.deleteDataRecord(atomID - 1);
	}

	/**
	 * Convenience method to set an Integer value in a data record.
	 * 
	 * @param field
	 * @param index
	 * @param value
	 */
	public void setInteger(String field, int index, int value)
	{
		setString(field, index, Integer.toString(value));
	}

	/**
	 * Shifts residue IDs of residues with specified names by some offset. You would
	 * want to do this if you merged multiple chains into one CT block and residues
	 * now have conflicting IDs. This method is meant to be called on a CT block with
	 * an m_atom child node. I guess this is an artificial restriction. Oh well.
	 *  
	 * @param byHowMuch Offset to add
	 * @param residueNames
	 * @throws MaeNodeException 
	 */
	public void shiftResidueIDs(int byHowMuch, String... residueNames) throws MaeNodeException
	{
		MaeNode m_atom = getNode("m_atom");
		if(m_atom == null)
			throw new MaeNodeException("shiftResidueIDs: This node doesn't have an m_atom block.");
		
		ArrayList<String> realResidueNames = new ArrayList<String>();
		for(String residueName : residueNames)
			realResidueNames.add(residueName);
		
		int numAtoms = m_atom.numDataRecords();
		for(int atomIndex = 0; atomIndex < numAtoms; atomIndex++)
		{
			String resName = m_atom.getString("s_m_pdb_residue_name", atomIndex).trim();
			if(realResidueNames.contains(resName))
			{
				int resNum = m_atom.getInteger("i_m_residue_number", atomIndex);
				m_atom.setInteger("i_m_residue_number", atomIndex, resNum + byHowMuch);
			}
		}
	}

	/**
	 * Copies PBC parameters from another CT blcok to this one.
	 * 
	 * @param from
	 * @param to
	 */
	public void copyPBCParameters(MaeNode from)
	{
		if(from.getColumn("r_chorus_box_ax") != null)
		{
			String[] pbcHeadings = { "r_chorus_box_ax", "r_chorus_box_ay", "r_chorus_box_az",
					"r_chorus_box_bx", "r_chorus_box_by", "r_chorus_box_bz", 
					"r_chorus_box_cx", "r_chorus_box_cy", "r_chorus_box_cz" };
			
			for(String key : pbcHeadings)
				this.setString(key, from.getString(key, 0));
		}
	}
	
	private String commaSeparatedList(String... items)
	{
		StringBuffer out = new StringBuffer();
		
		for(int i = 0; i < items.length - 1; i++)
			out.append(items[i] + ", ");
		
		out.append(items[items.length - 1]);		
		return out.toString();
	}

	public void dataRecordsAsSQLToStream(PrintStream out)
	{
		out.print("CREATE TABLE `" + getName() + "` (");
		
		ArrayList<String> foo = new ArrayList<String>();
		foo.add("`id` INTEGER PRIMARY KEY");
		
		for(String heading : headings_)
		{
			String type = "VARCHAR(255)";
			if(heading.startsWith("i_")) // integer
				type = "INTEGER";
			else if(heading.startsWith("r_")) // real number
				type = "DOUBLE";
			else if(heading.startsWith("s_")) // string
				type = "VARCHAR(255)";
			foo.add("`" + heading + "` " + type);
		}
		String[] bar = {};
		out.println("(" + commaSeparatedList(foo.toArray(bar)) + ");");

		
		out.print("INSERT INTO `" + getName() + "` (");
		//for(String heading : headings_)
		//	out.print("`");
		
		for (int recordID = 0; recordID < numDataRecords(); recordID++)
		{
			for(String heading : headings_)
			{
				
			}
		}
		
	}
	
}
