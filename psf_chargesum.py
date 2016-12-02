#!/usr/bin/env python
#
# Returns the overall charge sum of a PSF-format system
import argparse
import MDAnalysis as mda

def main():
    ap = argparse.ArgumentParser(description='Calculate sum of all charges in a PSF')
    ap.add_argument('psf')
    args = ap.parse_args()

    u = mda.Universe(args.psf)

    charge_sum = 0.0
    for atom in u.atoms:
        charge_sum += atom.charge
    
    print('Total charge for %s is %.f' % (args.psf, charge_sum))

if __name__ == '__main__':
    import sys
    sys.exit(int(main() or 0))