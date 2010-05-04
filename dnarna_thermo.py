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
import math

params_G = {
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
    'GrayINN2': { 'AA': -0.97,
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

# ∆H parameters
params_H = {
# Table I(B) - INN model with one nonunique initiation factor
    'GrayINN': {  'AA': -6.62,
                  'UU': -4.86,
                  'CC': -6.06,
                  'GG': -10.74,
                  'AU': -7.14,
                  'UA': -9.50,
                  'AC': -7.82,
                  'CA': -7.46,
                  'AG': -7.27,
                  'GA': -9.24,
                  'UG': -7.59,
                  'GU': -8.81,
                  'UC': -8.09,
                  'CU': -3.13,
                  'GC': -10.29,
                  'CG': -12.37,
                  'Initiation': -4.55 },
            
    # Table I(B) - Parameters from Sugimoto et al. with one nonunique initiation factor
    'Sugimoto': { 'AA': -7.8,
                  'UU': -11.5,
                  'CC': -9.3,
                  'GG': -12.8,
                  'AU': -8.3,
                  'UA': -7.8,
                  'AC': -5.9,
                  'CA': -9.0,
                  'AG': -9.1,
                  'GA': -5.5,
                  'UG': -10.4,
                  'GU': -7.8,
                  'UC': -8.6,
                  'CU': -7.0,
                  'GC': -8.0,
                  'CG': -16.3,
                  'Initiation': 1.9 },
                  
    # Table I(C) - Alternative INN parameters with end-factors
    'GrayINN2': { 'AA': -6.56,
                  'UU': -5.32,
                  'CC': -5.82,
                  'GG': -10.29,
                  'AU': -8.32,
                  'UA': -8.37,
                  'AC': -7.79,
                  'CA': -7.27,
                  'AG': -7.88,
                  'GA': -8.66,
                  'UG': -6.90,
                  'GU': -9.11,
                  'UC': -6.91,
                  'CU': -4.09,
                  'GC': -9.63,
                  'CG': -12.01,
                  'EA': -2.42, # E means end
                  'AE': -2.74,
                  'EU': -1.83,
                  'UE': -1.16,
                  'EG': -3.38,
                  'GE': -2.78,
                  'EC': -2.99,
                  'CE': -3.95 }
}

# ∆S parameters - caution! Units are cal/(mol*K)
params_S = {
# Table I(B) - INN model with one nonunique initiation factor
    'GrayINN': {  'AA': -18.3,
                  'UU': -14.6,
                  'CC': -13.4,
                  'GG': -26.1,
                  'AU': -20.7,
                  'UA': -28.2,
                  'AC': -19.5,
                  'CA': -20.0,
                  'AG': -18.5,
                  'GA': -23.9,
                  'UG': -20.1,
                  'GU': -24.0,
                  'UC': -21.9,
                  'CU': -6.7,
                  'GC': -24.7,
                  'CG': -34.1,
                  'Initiation': -24.0 },
            
    # Table I(B) - Parameters from Sugimoto et al. with one nonunique initiation factor
    'Sugimoto': { 'AA': -21.9,
                  'UU': -36.4,
                  'CC': -23.2,
                  'GG': -31.9,
                  'AU': -23.9,
                  'UA': -23.2,
                  'AC': -12.3,
                  'CA': -26.1,
                  'AG': -23.5,
                  'GA': -13.5,
                  'UG': -28.4,
                  'GU': -21.6,
                  'UC': -22.9,
                  'CU': -19.7,
                  'GC': -17.1,
                  'CG': -47.1,
                  'Initiation': -3.9 },
                  
    # Table I(C) - Alternative INN parameters with end-factors
    'GrayINN2': { 'AA': -18.0,
                  'UU': -16.2,
                  'CC': -12.4,
                  'GG': -24.5,
                  'AU': -24.5,
                  'UA': -24.5,
                  'AC': -20.0,
                  'CA': -18.9,
                  'AG': -20.0,
                  'GA': -22.2,
                  'UG': -17.6,
                  'GU': -25.3,
                  'UC': -18.8,
                  'CU': -9.3,
                  'GC': -23.4,
                  'CG': -31.9,
                  'EA': -12.8, # E means end
                  'AE': -13.9,
                  'EU': -10.0,
                  'UE': -8.2,
                  'EG': -15.7,
                  'GE': -14.3,
                  'EC': -14.9,
                  'CE': -17.0 }
}

def dnarna_annealing_val(params, seq):
    """Calculates thermodynamic value based on nearest-neighbor parameter set. """
    s = seq.upper()
    val = 0
    if 'Initiation' in params:
        val += params['Initiation']
    
    first = "E%s" % seq[0]
    last = "%sE" % seq[-1]
    if first in params:
        val += params[first]
    if last in params:
        val += params[last]
    
    for i in xrange(len(s) - 1):
        bases = s[i:i+2]
        try:
            energy = params[bases]
            # print "%s contributes %.2f" % (bases, energy)
            val += energy
        except KeyError:
            print "HEY! No parameter found for %s so don't trust the result" % bases
        
    return val


if __name__ == '__main__':
    seq = sys.argv[1]
    # T = (273.15+37)
    T = 300
    R = 1.9858775 / 1000.0
    conc = 80. / 1e6 # Concentration. Units are M, I guess
    print "For annealing of DNA/RNA hybrid with RNA sequence 5' %s 3':" % seq.upper()
    for param_set_name in params_G:
        dg = dnarna_annealing_val(params_G[param_set_name], seq)
        dh = dnarna_annealing_val(params_H[param_set_name], seq)
        # Units of ds are cal/(mol*K) so convert to kcal/(mol*K)
        ds = dnarna_annealing_val(params_S[param_set_name], seq) / 1000.0
        print "%s: ∆G = %.1f kcal/mol" % (param_set_name, dg)
        print "%s: ∆H = %.1f kcal/mol" % (param_set_name, dh)
        print "%s: ∆S = %.1f kcal/mol" % (param_set_name, ds)
        print "%s: At %.0fK, ∆H - T∆S = %.1f kcal/mol" % (param_set_name, T, dh - T*ds)
        # Tm = dh / (ds + R * math.log(conc))
        Tm = 1.0 / (R * math.log(conc/4) / dh + ds / dh)
        print "%s: For single strand concentration %.2fµM, Tm = %.1fK" % (param_set_name, conc*1e6, Tm)
        