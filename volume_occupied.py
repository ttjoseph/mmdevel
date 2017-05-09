#!/usr/bin/env python
import mdtraj as md
import numpy as np
import argparse
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Returns the volume occupied by a ligand over the course of the simulation.')
    parser.add_argument('ligand_resname', help='Residue name of ligand you are interested in (e.g. SKE)')
    parser.add_argument('psf', help='Topology file, in PSF format')
    parser.add_argument('dcd', nargs='+', help='Trajectory files')
    args = parser.parse_args()

    # Load the trajectory, superpose according to the protein in the first frame
    print >>sys.stderr, 'Loading and superposing trajectory: %s' % args.dcd
    traj = md.load(*args.dcd, top=args.psf)
    traj.superpose(traj, atom_indices=traj.topology.select('protein'))

    # Determine the max and min coordinates of the ligand over the trajectory
    print >>sys.stderr, 'Processing frames...'
    atoms = traj.topology.select('resname %s' % args.ligand_resname)
    if len(atoms) == 0:
        print >>sys.stderr, 'Cannot find resname %s in the topology you provided. Is that really the ligand you want?' % args.ligand_resname
        exit(1)
    max_x, max_y, max_z = -10000, -10000, -10000
    min_x, min_y, min_z = 10000, 10000, 10000
    for frame in range(traj.n_frames):
        for atom in atoms:
            x, y, z = traj.xyz[frame, atom, :]
            max_x = x if x > max_x else max_x
            max_y = y if y > max_y else max_y
            max_z = z if z > max_z else max_z
            min_x = x if x < min_x else min_x
            min_y = y if y < min_y else min_y
            min_z = z if z < min_z else min_z

    # Calculate the radius of the smallest sphere that could contain the ligand during the trajectory
    x_width = max_x - min_x
    y_width = max_y - min_y
    z_width = max_z - min_z

    sphere_radius = np.max(np.array([x_width, y_width, z_width])) / 2
    # Units are nanometers internally to MDTraj
    print 'Box widths (A): %f %f %f' % (x_width * 10, y_width * 10, z_width * 10)
    print 'Sphere radius (A): %.2f' % (sphere_radius * 10)

