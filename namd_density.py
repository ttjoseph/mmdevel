#!/usr/bin/env python
import argparse
import MDAnalysis as mda

CM3_PER_ANGSTROM3 = 1.0e-24
G_PER_AMU = 1.66054e-24

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Calculate mass density (in g/cm^3) over time from NAMD PSF and log')
    ap.add_argument('psf', help='NAMD PSF file, which we will use to get atom masses')
    ap.add_argument('log', nargs='+', help='NAMD output log file(s) (not the trajectories)')
    args = ap.parse_args()

    # Calculate total mass of all atoms in system
    u = mda.Universe(args.psf)
    total_mass = 0.0
    for a in u.atoms:
        total_mass += a.mass

    total_mass *= G_PER_AMU

    ts_col, volume_col = -1, -1
    for log in args.log:
        f = open(log, 'r')
        for line in f:
            # I don't know which column to look for. Let the computer figure it out
            if volume_col < 0 and line.startswith('ETITLE:'):
                tokens = line.split()
                ts_col = tokens.index('TS')
                volume_col = tokens.index('VOLUME')

            if volume_col >= 0 and line.startswith('ENERGY:'):
                tokens = line.split()
                ts = int(tokens[ts_col])
                volume = float(tokens[volume_col]) * CM3_PER_ANGSTROM3
                print ts, total_mass / volume

        f.close()
