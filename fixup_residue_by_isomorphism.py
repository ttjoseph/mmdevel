#!/usr/bin/env python
#
# Given a known-good residue template and the same residue but with weird atom names
# (perhaps from an experimentalist), update the known-good residue atom coordinates
# from the corresponding atoms in the weirdly-named residue.
#
# Only supports the first residue of each PDB, so you probably want to just work with
# single-residue PDBs.
#
# This allows you to easily convert, say, experimentally determined lipid configurations
# into simulation systems. Currently you have to do this one residue at a time, because
# I can't be bothered to make it, say, parse the entirety of CHARMM force field topologies
# and map all the residues.
import sys
import argparse
import MDAnalysis as mda
import networkx as nx

def residue_to_graph(u, residx, include_hydrogens=False):
    """Converts a residue in an MDAnalysis Universe, at 0-based index residx, to a NetworkX graph.
    
    Each node in the graph is identifed by a tuple of atom name and atom index. Edges come from
    the bond information included in the Universe.

    Hydrogens are ignored by default, unless you specify include_hydrogens=True.
    """
    # MDAnalysis does not make the guessed bonds easy to retrieve.
    # For example, you can't get the bonds out of a particular residue - it says
    # there is not bond information in the Universe even if there actually is.
    # So, we have to make a list of atoms in the residue of interest, then
    # manually select out the bonds between only those atoms, so we can make
    # our graph.
    # Perhaps there is a better way but the MDAnalysis documentation sure doesn't
    # make it easy to find out.

    # Make it optional to ignore hydrogens, since they may not be present in the experimental structure
    # include_hydrogens is False and name starts with H:          not (1 and 1) -> 0 = don't use
    # include_hydrogens is False and name does not start with H:  not (1 and 0) -> 1 = use
    # include_hydrogens is True and name starts with H:           not (0 and 1) -> 1 = use
    # include_hydrogens is True and name does not start with H:   not (0 and 0) -> 1 = use
    def should_use(a):
        return not (include_hydrogens is False and a.name.startswith('H'))

    # Get a list of atom indices
    atom_ixs = [a.ix for a in u.residues[residx].atoms if should_use(a)]
    atom_data = [(a.name, a.ix) for a in u.residues[residx].atoms if should_use(a)]
    # Pick out only the bonds for this residue
    bond_tuples = []
    for i in range(len(u.bonds)):
        a1, a2 = u.bonds[i].atoms
        if a1.ix in atom_ixs and a2.ix in atom_ixs and should_use(a1) and should_use(a2):
            bond_tuples.append(((a1.name, a1.ix), (a2.name, a2.ix)))
    g = nx.Graph()
    g.add_nodes_from(atom_data)
    g.add_edges_from(bond_tuples)    
    return g


def main():
    ap = argparse.ArgumentParser(description='Maps one residue with screwed up atom names to a gold standard residue, ignoring hydrogens')
    ap.add_argument('pdb_from', help='PDB of residue whose atom names are wrong')
    ap.add_argument('pdb_to', help='PDB of residue whose atom names are correct')
    ap.add_argument('pdb_out', help='Filename to write updated version of pdb_to')
    ap.add_argument('--include-hydrogens', '-i', action='store_true', help='Include hydrogens in mapping')
    ap.add_argument('--verbose', '-v', action='store_true', help='Be noisy and show the mapping')
    args = ap.parse_args()

    u_from = mda.Universe(args.pdb_from, guess_bonds=True)
    u_to = mda.Universe(args.pdb_to, guess_bonds=True)

    # Build graphs for each residue
    g_from = residue_to_graph(u_from, 0, args.include_hydrogens)
    g_to = residue_to_graph(u_to, 0, args.include_hydrogens)

    # Get the isomorphism between the two atom connectivities, complain if there isn't one
    matcher = nx.algorithms.isomorphism.GraphMatcher(g_from, g_to)
    if matcher.is_isomorphic() is False:
        print('Error: The first residue of each of those two PDBs are not isomorphic.', file=sys.stderr)
        exit(1)

    if args.verbose:
        for k, v in matcher.mapping.items():
            print(k, v, file=sys.stderr)

    u_from.residues[0].resname = u_to.residues[0].resname
    for atom_from, atom_to in matcher.mapping.items():
        a1_name, a1_idx = atom_from
        a2_name, a2_idx = atom_to
        # Check that the first character of each atom name is the same, as this customarily
        # indicates what element it is.
        if a1_name[0] != a2_name[0]:
            print(f'Warning: {a1_name} in {args.pdb_from} may not be the same element as {a2_name} in {args.pdb_to}', file=sys.stderr)
        u_to.residues[0].atoms[a2_idx].position = u_from.residues[0].atoms[a1_idx].position

    u_to.select_atoms('all').write(args.pdb_out)
    print(f'An MDAnalysis warning above about unit cell dimensions is probably OK.', file=sys.stderr)
    print(f'Wrote {args.pdb_out} with coordinates from the first residue in {args.pdb_from} and names from the first residue in {args.pdb_to}.', file=sys.stderr)
    print(f'If there were hydrogens in {args.pdb_to}, their coordinates were not modified!', file=sys.stderr)
    print(f'I suppose you will be doing some copy and paste next.', file=sys.stderr)


if __name__ == "__main__":
    main()
