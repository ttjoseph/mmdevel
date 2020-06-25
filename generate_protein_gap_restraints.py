#!/usr/bin/python
#
# Finds gaps in PDB files, which you probably want to fill in with a loop modeler or something
import argparse
import MDAnalysis as mda
from scipy.spatial.distance import euclidean

ap = argparse.ArgumentParser(description='Find missing loops in a protein structure')
ap.add_argument('pdb', help='Protein structure')
args = ap.parse_args()

u = mda.Universe(args.pdb)
protein = u.select_atoms('protein')

last_res = None

for res in protein.residues:
    if last_res is None:
        last_res = res

    # Mod 10000 because resid can only be four decimal digits
    if res.resid > last_res.resid + 1:
        print ('# Discontinuity between residues %s:%s%d and %s:%s%d (length %d)' \
               % (last_res.segid, last_res.resname, last_res.resid % 10000,
               res.segid, res.resname, res.resid % 10000, res.resid - last_res.resid - 1))
        colvar_name = 'gap_%s%d_%s%d' % (last_res.segid, last_res.resid % 10000, res.segid, res.resid % 10000)
        a1 = res.atoms.select_atoms('name CA')[0]
        a2 = last_res.atoms.select_atoms('name CA')[0]
        distance = euclidean(a1.position, a2.position)
        print('# Distance between CA: %f' % distance)

        print("""colvar {
    name %s
    distance {
        group1 { atomNumbers %d }
        group2 { atomNumbers %d }
    }
}""" % (colvar_name, a1.id, a2.id))
        print("""harmonic {
    colvars %s
    centers %f
    forceConstant 10.0
}""" % (colvar_name, distance))
        
    last_res = res
