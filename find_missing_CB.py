from pdb_to_psf import *

m = Molecule(sys.argv[1])
print >>sys.stderr, "This is a list of Gly residues (that don't have CB atoms)."
my_atoms = [a for a in m.atoms if a.resname == "GLY" and a.atomname == "N"]
first_resid = my_atoms[0].resid
for a in my_atoms:
	print (a.resid - first_resid + 1)
