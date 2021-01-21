#!/usr/bin/env python
#
# Given a known-good residue template and the same residue but with weird atom names
# (perhaps from an experimentalist), update the known-good residue heavy atom coordinates
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
from typing import Mapping
import MDAnalysis as mda
import networkx as nx
from networkx.algorithms.isomorphism.tree_isomorphism import tree_isomorphism

def residue_to_graph(u, residx):
    """Converts a residue in an MDAnalysis Universe, at 0-based index residx, to a NetworkX graph.
    
    Each node in the graph is identifed by a tuple of atom name and atom index. Edges come from
    the bond information included in the Universe.
    """
    # MDAnalysis does not make the guessed bonds easy to retrieve.
    # For example, you can't get the bonds out of a particular residue - it says
    # there is not bond information in the Universe even if there actually is.
    # So, we have to make a list of atoms in the residue of interest, then
    # manually select out the bonds between only those atoms, so we can make
    # our graph.
    # Perhaps there is a better way but the MDAnalysis documentation sure doesn't
    # make it easy to find out.

    # Get a list of atom indices
    atom_ixs = [a.ix for a in u.residues[residx].atoms if not a.name.startswith('H')]
    atom_data = [(a.name, a.ix) for a in u.residues[residx].atoms if not a.name.startswith('H')]
    # Pick out only the bonds for this residue
    bond_tuples = []
    for i in range(len(u.bonds)):
        a1, a2 = u.bonds[i].atoms
        if a1.ix in atom_ixs and a2.ix in atom_ixs and not a1.name.startswith('H') and not a2.name.startswith('H'):
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
    args = ap.parse_args()

    u_from = mda.Universe(args.pdb_from, guess_bonds=True)
    u_to = mda.Universe(args.pdb_to, guess_bonds=True)

    # Build graphs for each residue
    g_from = residue_to_graph(u_from, 0)
    g_to = residue_to_graph(u_to, 0)
    
    # Get the isomorphism between the two atom connectivities, complain if there isn't one
    mapping = tree_isomorphism(g_from, g_to)
    if len(mapping) == 0:
        print('Error: The first residue of each of those two PDBs are not isomorphic.', file=sys.stderr)
        exit(1)

    u_from.residues[0].resname = u_to.residues[0].resname
    for atom_from, atom_to in mapping:
        a1_name, a1_idx = atom_from
        a2_name, a2_idx = atom_to
        if a1_name[0] != a2_name[0]:
            print(f'Warning: {a1_name} in {pdb_from} may not be the same element as {a2_name} in {pdb_to}', file=sys.stderr)
        u_to.residues[0].atoms[a2_idx].position = u_from.residues[0].atoms[a1_idx].position

    u_to.select_atoms('resid 1').write(args.pdb_out)
    print(f'Wrote a single residue to {args.pdb_out}. I suppose you will be doing some copy and paste next.', file=sys.stderr)


if __name__ == "__main__":
    main()