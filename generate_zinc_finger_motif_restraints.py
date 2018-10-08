#!/usr/bin/python
#
# Generates restraints for zinc finger motifs
import argparse
import MDAnalysis as mda
from scipy.spatial.distance import euclidean

def emit_harmonic_restraint(colvar_name, a1_id, a2_id, distance, force_constant):
    print("""colvar {
    name %s
    distance {
        group1 { atomNumbers %d }
        group2 { atomNumbers %d }
    }
}""" % (colvar_name, a1_id, a2_id))
    print("""harmonic {
    colvars %s
    centers %f
    forceConstant %f
}""" % (colvar_name, distance, force_constant))


ap = argparse.ArgumentParser(description='Generate restraints for zinc finger motifs')
ap.add_argument('pdb', help='Protein coordinates')
ap.add_argument('--force-constant', default=10.0, help='Restraint force constant')
args = ap.parse_args()

u = mda.Universe(args.pdb)

# Find all zinc atoms
zincs = u.select_atoms('resname ZN2')

for zinc in zincs:
    zinc_seltext = 'atom %s %d %s' % (zinc.segid, zinc.resnum, zinc.name)
    # Find all cysteine SG atoms
    cys_sg = u.select_atoms('resname CYS and name SG and around 4.5 %s' % zinc_seltext)
    for sg in cys_sg:
        # Generate restraint
        distance = euclidean(sg.position, zinc.position)
        print('# Restraint for ZN2 (%s:%d) vs CYS (%s:%d.SG)' % (zinc.segid, zinc.resid, sg.segid, sg.resid))
        emit_harmonic_restraint('zn_%s%d_%s%d' % (zinc.segid, zinc.resid, sg.segid, sg.resid),
            zinc.id, sg.id, distance, args.force_constant)
        

    his_ne2 = u.select_atoms('resname HSD and name NE2 and around 4.5 %s' % zinc_seltext)
    for ne2 in his_ne2:
        distance = euclidean(ne2.position, zinc.position)
        print('# Restraint for ZN2 (%s:%d) vs HSD (%s:%d.NE2)' % (zinc.segid, zinc.resid, ne2.segid, ne2.resid))
        emit_harmonic_restraint('zn_%s%d_%s%d' % (zinc.segid, zinc.resid, ne2.segid, ne2.resid),
            zinc.id, ne2.id, distance, args.force_constant)
