#!/usr/bin/env python
import argparse
import sys
import numpy as np
import MDAnalysis as mda

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
    ap.add_argument('helixlist', help='Helix residue ID ranges, inclusive, separated by commas with no whitespace; e.g. 10-20,60-90' )
    ap.add_argument('psf', help='PSF topology file')
    ap.add_argument('dcd', nargs='+', help='DCD trajectory file[s]')
    args = ap.parse_args()

    u = mda.Universe(args.psf, args.dcd)

    # Parse the helix residue ID ranges
    helices, labels = [], []
    for helix in args.helixlist.split(','):
        resids = [int(x) for x in helix.split('-')]
        if len(resids) != 2:
            raise ValueError('Your helix specification %s does not have exactly two residue IDs' % helix)
        helices.append(u.select_atoms('protein and name CA and resid %d:%d' % (resids[0], resids[1])))
        labels.append('"Radial %d-%d"' % (resids[0], resids[1]))
        labels.append('"Lateral %d-%d"' % (resids[0], resids[1]))
    # Print CSV label
    print ','.join(labels)

    protein = u.select_atoms('protein and name CA')
    for ts in u.trajectory:
        c1, c2, c3 = protein.principal_axes()
        # Force the first axis to always point up
        if c1[2] < 0: c1 = -c1

        # Calculate angles between the various principal axes

        # The vector in the xy plane originating at the channel COM and extending to the COM of the TM2 vector
        # is normal to the plane against which we measure the TM2 tilting angle for radial tilt.
        vals = []
        for helix in helices:
            a1, a2, a3 = helix.principal_axes()
            # Force the first axis to point up, always
            if a1[2] < 0: a1 = -a1
            radial_plane_normal = helix.center_of_mass() - protein.center_of_mass()
            radial_plane_normal[2] = 0
            radial_tilting_angle = 90-angle_in_degrees(a1, radial_plane_normal)

            # The plane relative to which we measure lateral tilting angle is orthogonal to the radial angle plane
            lateral_plane_normal = np.cross(radial_plane_normal, c1)
            lateral_tilting_angle = 90-angle_in_degrees(a1, lateral_plane_normal)

            vals.append(radial_tilting_angle)
            vals.append(lateral_tilting_angle)

            # Print out diagnostic VMD arrow drawing commands if we see a really wacky lateral tilting angle
            # if lateral_tilting_angle > 20:
            #     logging.info("Frame %d, radial=%.2f, lateral=%.2f" % (frame_num, radial_tilting_angle,
            #                                                           lateral_tilting_angle))
            #     print "graphics top color blue"
            #     vmd_draw_arrow(channel.centerOfMass(), c1*20)
            #     vmd_draw_arrow(channel.centerOfMass(), radial_plane_normal)
            #     print "graphics top color red"
            #     vmd_draw_arrow(tm.centerOfMass(), a1*20)
            #     exit()

        print ','.join(['%.2f' % x for x in vals])

