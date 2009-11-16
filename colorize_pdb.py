from AMBER import *
import sys

info("Hi there.")

if len(sys.argv) < 3:
    print >>sys.stderr, "Usage: <pdb> <value-filename>"
    sys.exit()

pdb = PDB(sys.argv[1])
pdb.nuke_solvent()


for line in open(sys.argv[2]):
    try:
        values = [float(x) for x in line.split()]
    except ValueError:
        print "That be wack, yo."
        sys.exit()

    # Color the PDB by residue using occupancy and tempfactor fields
    for atom in pdb.atoms:
        val = values[atom.resid - 1]
        if val > 5: val = 5
        if val < -1: val = -5
        atom.occupancy = atom.tempfactor = val

    pdb.write(sys.stdout)