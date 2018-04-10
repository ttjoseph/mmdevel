#!/usr/bin/env python
import argparse
import sys
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt

# Convenient resid ranges for the transmembrane helices of selected proteins
# 5TVN: 5HT2B receptor
# 5C1M: MuOR
# 4DJH: KappaOR
HELIX_LISTS = {'5tvn': '50-80,90-117,127-159,169-193,215-245,321-350,356-382',
               '5c1m': '65-98,102-130,137-169,181-205,226-259,270-304,313-343',
               '4djh': '56-86,93-121,128-161,171-196,219-255,270-299,310-332'}


def closest_to_z_axis(*vectors):
    """Of the vectors supplied, returns the one that has the smallest angle with respect to the Z axis."""
    z_axis = np.array([0, 0, 1])
    min_angle = 900
    min_vector = None
    for v in vectors:
        angle = angle_in_degrees(v, z_axis)
        # TODO: This doesn't work because we should consider positive and negative vectors that are parallel the same
        if np.abs(angle) < min_angle:
            min_vector = v
            min_angle = angle
    return min_vector


def angle_in_degrees(u, v):
    """Returns the angle between two vectors, in degrees."""
    angle = np.arccos(np.dot(u, v)/np.linalg.norm(u)/np.linalg.norm(v))
    return angle*(180.0/np.pi)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Calculate helix tilts over a trajectory, relative to longest protein axis')
    ap.add_argument('helixlist', help='Helix residue ID ranges, inclusive, separated by commas with no whitespace; e.g. 10-20,60-90. Or a list alias.' )
    ap.add_argument('psf', help='PSF topology file')
    ap.add_argument('dcd', nargs='+', help='DCD trajectory file[s]')
    ap.add_argument('--absolute', action='store_true', help='Plot absolute angles rather than those relative to first frame')
    args = ap.parse_args()

    u = mda.Universe(args.psf, args.dcd)
    out_file = open('helix-tilt.csv', 'w')

    # Parse the helix residue ID ranges
    helixlist = HELIX_LISTS[args.helixlist] if args.helixlist in HELIX_LISTS else args.helixlist
    helices, labels = [], []
    first_radial_tilting_angle, first_lateral_tilting_angle = {}, {}
    for helix in helixlist.split(','):
        resids = [int(x) for x in helix.split('-')]
        if len(resids) != 2:
            raise ValueError('Your helix specification %s does not have exactly two residue IDs' % helix)
        helices.append(u.select_atoms('protein and name CA and resid %d:%d' % (resids[0], resids[1])))
        labels.append('"Rad%d-%d"' % (resids[0], resids[1]))
        labels.append('"Lat%d-%d"' % (resids[0], resids[1]))
        first_radial_tilting_angle[helices[-1]] = None
        first_lateral_tilting_angle[helices[-1]] = None
    # Print CSV label
    print >>out_file, ','.join(labels)

    # Iterate over frames in the trajectory
    protein = u.select_atoms('protein and name CA')
    for ts in u.trajectory:
        c3, c2, c1 = protein.principal_axes()
        # Force the major axis to always point up. Normally this should be the z-axis
        # print >>sys.stderr, c1
        if c1[2] < 0: c1 = -c1
        # Assert that we've chosen the pricipal axis closest to the Z axis
        assert(c1[2] > c1[0])
        assert(c1[2] > c1[1])

        # Calculate angles between the various principal axes.
        # The vector in the xy plane originating at the channel COM and extending to the COM of the TM2 vector
        # is normal to the plane against which we measure the TM2 tilting angle for radial tilt.
        vals = []
        for helix in helices:
            a3, a2, a1 = helix.principal_axes()
            # Force the first axis to point up, always
            if a1[2] < 0: a1 = -a1
            # Assert that we've chosen the pricipal axis closest to the Z axis
            assert (a1[2] > a1[0])
            assert (a1[2] > a1[1])
            radial_plane_normal = helix.center_of_mass() - protein.center_of_mass()
            # print >>sys.stderr, radial_plane_normal[2]
            # radial_plane_normal[2] = 0
            radial_tilting_angle = 90-angle_in_degrees(a1, radial_plane_normal)
            if first_radial_tilting_angle[helix] is None:
                first_radial_tilting_angle[helix] = radial_tilting_angle

            # The plane relative to which we measure lateral tilting angle is orthogonal to the radial angle plane
            lateral_plane_normal = np.cross(radial_plane_normal, c1)
            lateral_tilting_angle = 90-angle_in_degrees(a1, lateral_plane_normal)
            if first_lateral_tilting_angle[helix] is None:
                first_lateral_tilting_angle[helix] = lateral_tilting_angle

            # If specified in command line flag, keep and subtract the angles from the first frame
            if args.absolute is False:
                radial_tilting_angle -= first_radial_tilting_angle[helix]
                lateral_tilting_angle -= first_lateral_tilting_angle[helix]

            vals.append(radial_tilting_angle)
            vals.append(lateral_tilting_angle)

        print >>out_file, ','.join(['%.2f' % x for x in vals])


    # plt.figure()
    # plt.savefig('helices.png', bbox_inches='tight')
    print """# Gnuplot script to plot this data
set datafile separator ','
set xlabel 'Radial tilt (degrees)'
set ylabel 'Lateral tilt (degrees)'
set term png
fname = 'helix-tilt.csv'
set output 'helix-tilt.png'
"""
    s = []
    for i in range(0, len(labels), 2):
        # Using two columns at a time
        s.append("fname using %d:%d title 'TM%d'" % (i+1, i+2, i/2+1))

    print 'plot %s' % ', '.join(s)
