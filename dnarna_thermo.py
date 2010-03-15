#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Estimates delta-G for a DNA:RNA oligo given the RNA sequence in 5' to 3' direction.
# Uses parameters from:
#
# Gray DM. Derivation of Nearest Neighbor Properties from Data on Nucleic
# Acid Oligomers. II. Thermodynamic Parameters of DNA:RNA Hybrids and DNA Duplexes.
# Biopolymers, Vol 42, 795-810 (1997).
#
# Sugimoto N et al. Biochemistry 34, 11211-11216 (1995).
#
# Tom Joseph <thomas.joseph@mssm.edu>
import sys

params = {
# Table I(B) - INN model with one nonunique initiation factor
    'GrayINN': {  'AA': -0.95,
                  'UU': -0.34,
                  'CC': -1.92,
                  'GG': -2.64,
                  'AU': -0.71,
                  'UA': -0.75,
                  'AC': -1.76,
                  'CA': -1.25,
                  'AG': -1.54,
                  'GA': -1.82,
                  'UG': -1.37,
                  'GU': -1.36,
                  'UC': -1.28,
                  'CU': -1.07,
                  'GC': -2.62,
                  'CG': -1.78,
                  'Initiation': 2.90 },
            
    # Table I(B) - Parameters from Sugimoto et al. with one nonunique initiation factor
    'Sugimoto': { 'AA': -1.0,
                  'UU': -0.2,
                  'CC': -2.1,
                  'GG': -2.9,
                  'AU': -0.9,
                  'UA': -0.6,
                  'AC': -2.1,
                  'CA': -0.9,
                  'AG': -1.8,
                  'GA': -1.3,
                  'UG': -1.6,
                  'GU': -1.1,
                  'UC': -1.5,
                  'CU': -0.9,
                  'GC': -2.7,
                  'CG': -1.7,
                  'Initiation': 3.1 },
                  
    # Table I(C) - Alternative INN parameters with end-factors
    'GrayINN2': { 'AA': -0.7,
                  'UU': -0.30,
                  'CC': -1.98,
                  'GG': -2.69,
                  'AU': -0.73,
                  'UA': -0.77,
                  'AC': -1.57,
                  'CA': -1.42,
                  'AG': -1.68,
                  'GA': -1.79,
                  'UG': -1.45,
                  'GU': -1.26,
                  'UC': -1.10,
                  'CU': -1.21,
                  'GC': -2.36,
                  'CG': -2.12,
                  'EA': 1.55, # E means end
                  'AE': 1.56,
                  'EU': 1.28,
                  'UE': 1.39,
                  'EG': 1.49,
                  'GE': 1.66,
                  'EC': 1.62,
                  'CE': 1.34 }
}

def dnarna_annealing_deltag(params, seq):
    """Calculates ∆G """
    s = seq.upper()
    dg = 0
    if 'Initiation' in params:
        dg += params['Initiation']
    
    first = "E%s" % seq[0]
    last = "%sE" % seq[-1]
    if first in params:
        dg += params[first]
    if last in params:
        dg += params[last]
    
    for i in xrange(len(s) - 1):
        bases = s[i:i+2]
        try:
            energy = params[bases]
            # print "%s contributes %f" % (bases, energy)
            dg += energy
        except KeyError:
            print "HEY! No parameter found for %s so don't trust the result" % bases
        
    return dg


if __name__ == '__main__':
    seq = sys.argv[1]
    print "For sequence %s:" % seq.upper()
    for name in params:
        dg = dnarna_annealing_deltag(params[name], seq)
        print "%s ∆G = %.1f kcal/mol" % (name, dg)