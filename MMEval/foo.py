import AMBER
import mmeval

m = AMBER.AmberSystem("complex.top.x", "complex.crd.1627")
num_atom_types = m.blocks['POINTERS'][1]
print "Number of atom types: %d" % num_atom_types
print "NONBONDED_PARM_INDEX length: %d" % len(m.blocks['NONBONDED_PARM_INDEX'])

energy = mmeval.bondenergy(m)
print "Ebond = %f" % energy
energy = mmeval.elecenergy(m)
print "Eelec = %f" % energy
energy = mmeval.vdwenergy(m)
print "Evdw = %f" % energy
