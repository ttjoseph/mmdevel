from pdb_to_psf import *

m = Molecule(sys.argv[1])
print >>sys.stderr, "You want %s atoms." % sys.argv[2]
my_atoms = [a for a in m.atoms if a.atomname == sys.argv[2]]
first_resid = my_atoms[0].resid
for a in my_atoms:
    print "%d %f %f %f" % (a.resid - first_resid + 1, a.x, a.y, a.z)
