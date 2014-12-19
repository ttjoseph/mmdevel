#!/usr/bin/env python 
import logging
import numpy 
from MDAnalysis import *
from MDAnalysis.analysis.distances import *

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

DISTANCE_CUTOFF = 15.0 # Angstroms

RESIDUES_PER_SUBUNIT = 311
# I202, V242, T255 as described in Nury et al.
IMPORTANT_RESIDUES_1 = (202, 242, 255)
IMPORTANT_RESIDUES_1a = (305+4, -1)
IMPORTANT_RESIDUES_2 = (209+4, 292+4, 265+4)

Structure = 'withxe_md3.gro'
Trajectory = ['withxe_md3.xtc', 'withxe_md4.xtc']

logging.info("Loading structure %s and trajectory %s..." % (Structure, Trajectory))
u = Universe(Structure, Trajectory)

logging.info("There are %d frames in the trajectory." % u.trajectory.numframes)

backbone = u.selectAtoms('name BB')
ligand = u.selectAtoms('name XE')

# Iterate thorugh each frame
pixels_out = 0
counts = numpy.zeros(RESIDUES_PER_SUBUNIT)
for ts in u.trajectory:
    # For each backbone atom, we want to find the minimum distance to a ligand
    d = distance_array(backbone.coordinates(), ligand.coordinates())
    # Each element of d is an array of distances to all ligands
    resid = 1
    for bb in d:
        val = min(bb)
        # We add 4 because the PDB 3P50 starts at resid 5. Thus our resid 1 is the real resid 5
        if val < DISTANCE_CUTOFF:
            counts[(resid-1)%RESIDUES_PER_SUBUNIT] += 1
        resid += 1

for i in xrange(0, len(counts)):
    # The GLIC subnits each start at resid 5
    # print "%d %d" % (i + 1 + 4, counts[i])
    print "%d " % (counts[i]),
