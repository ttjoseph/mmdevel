#!/usr/bin/env python
#
# Replaces MARTINI single-point waters ('W') with polarizable MARTINI waters ('PW').
# This involves adding the WP and WM sites. Only outputs the waters, which you will have
# to manually cut and paste into a PDB also containing your solute, and then put that
# through grompp to produce your TPR file.

import sys
import logging
import MDAnalysis
import numpy as np

class AtomRecord:
    '''Represents a single atom.'''
    def __init__(self, atom):
        self.atomid = atom.number
        self.atomname = atom.name
        self.resname = atom.residue.name
        self.chain = atom.segid
        self.resid = atom.resid
        self.x = atom.pos[0]
        self.y = atom.pos[1]
        self.z = atom.pos[2]
        self.occupancy = 0.0
        self.tempfactor = 0.0

    def write(self, fp):
        """Writes this atom, in PDB format."""
        
        fp.write("ATOM  %5d %4s %3s %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n" %
                 (self.atomid, self.atomname, self.resname, self.resid,
                 self.x, self.y, self.z, self.occupancy, self.tempfactor))

if __name__ == '__main__':
    logging.basicConfig(format='# %(message)s', level=logging.INFO)
    logging.info("Converts regular MARTINI waters to polarizable waters - Tom Joseph <thomas.joseph@mountsinai.org>")

    if len(sys.argv) < 2:
        logging.error("Usage: <gro-file>\n\nSends resulting PDB to stdout.")
        exit(1)

    u = MDAnalysis.Universe(sys.argv[1])
    waters = u.selectAtoms('resname W')

    logging.info("OK, I'm going to print only a water block.")

    count = 1
    for wat in waters:
        # We want to make a new PW residue, with a W atom at the center and
        # also a WP and WM atom, each of which are 0.14 nm from the W atom
        w = AtomRecord(wat)
        wp = AtomRecord(wat)
        wm = AtomRecord(wat)

        wp.x += 1.4  # Angstroms, not nm
        wm.y += 1.4

        w.resname = wp.resname = wm.resname = "PW"
        w.atomname = "W"
        wp.atomname = "WP"
        wm.atomname = "WM"

        w.atomid = count
        wp.atomid = count+1
        wm.atomid = count+2

        w.write(sys.stdout)
        wp.write(sys.stdout)
        wm.write(sys.stdout)

        count += 3
