#!/usr/bin/env python
#
# Shamelessly ripped off from:
#     http://www.ks.uiuc.edu/Research/vmd/script_library/scripts/orient/
import logging

from scipy import linalg
import MDAnalysis
import numpy as np

def inertia_tensor(atoms):
    """Computes the inertia tensor, whatever that is. Does not do any mass weighting."""
    # Remove center of mass
    coords = atoms.coordinates() - atoms.centerOfMass()
    ixx = ixy = ixz = iyy = iyz = izz = 0.0

    for c in coords:
        x, y, z = c[0], c[1], c[2]
        ixx += y*y+z*z
        ixy -= x*y
        ixz -= x*z
        iyy += x*x+z*z
        iyz -= y*z
        izz += x*x+y*y

    return np.array([[ixx, ixy, ixz], [ixy, iyy, iyz], [ixz, iyz, izz]])


def principal_axes(atoms):
    """Returns a matrix for which the columns are the principal axes of the specified subset of atoms.

    By definition there are three principal axes for a 3D rigid body. Vectors are returned in order of increasing
    eigenvalue, and apparently the eigenvector with the smallest eigenvalue is the long axis. So this works best
    with subsets of atoms that have one clear long axis, such as alpha helices."""
    # We can guarantee the inertia tensor is real and symmetric, so we can use linalg.eigh
    w, v = linalg.eigh(inertia_tensor(atoms))
    # Vectors are returned in order of increasing eigenvalue
    return v[:, 0], v[:, 1], v[:, 2]


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


def vmd_draw_arrow(o, v):
    print "draw arrow {%f %f %f} {%f %f %f}" % (o[0], o[1], o[2], o[0]+v[0], o[1]+v[1], o[2]+v[2])


if __name__ == '__main__':
    logging.basicConfig(format='# %(message)s', level=logging.INFO)
    logging.info("Computes helix angles and whatever - Tom Joseph <thomas.joseph@mountsinai.org>")

    # BasePath = '/Users/tom/Desktop/GLIC-CG/WithPFL/Intra-1/Restrained/'
    # orig = MDAnalysis.Universe(BasePath+'min_pfl.gro')
    # u = MDAnalysis.Universe(BasePath+'pr2.gro', BasePath+'md1-200kJ.xtc')

    BasePath = '/Users/tom/Desktop/GLIC-CG/EN_Sandbox2/'
    orig = MDAnalysis.Universe(BasePath+'min.gro')
    u = MDAnalysis.Universe(BasePath+'md1.gro', BasePath+'md1.xtc')

    # Atom indices below are valid only for the CG system.
    # TM2 A: bynum 510:562
    orig_tm2 = orig.selectAtoms('name BB and bynum 510:562')
    o1, o2, o3 = principal_axes(orig_tm2)
    # We ignore the extracellular domain when defining the channel because it flails all over the place
    channel = u.selectAtoms('name BB and resid 221:245')
    whole_channel = u.selectAtoms('name BB and resid 1:311')
    logging.info("There are %d backbone atoms in the channel." % len(whole_channel))
    if len(whole_channel) != 1555:
        logging.fatal("But there should be 1555 atoms. Something is wrong!")
        exit(1)
    bilayer = u.selectAtoms('resname POPC')
    # There are five TM2 helices, one for each subunit
    # TM2 is resid 221:245
    tm2 = []
    tm2.append(u.selectAtoms('name BB and bynum 510:562'))
    tm2.append(u.selectAtoms('name BB and bynum 1226:1278'))
    tm2.append(u.selectAtoms('name BB and bynum 1942:1994'))
    tm2.append(u.selectAtoms('name BB and bynum 2658:2710'))
    tm2.append(u.selectAtoms('name BB and bynum 3374:3426'))

    logging.info("There are %d frames in %s." % (u.trajectory.numframes, u.trajectory.filename))

    frame_num = 0
    z_axis = np.array([0, 0, 1])
    for ts in u.trajectory:
        frame_num += 1
        # if frame_num < 6770: continue
        c1, c2, c3 = principal_axes(channel)
        # Force the first axis to always point up
        if c1[2] < 0: c1 = -c1
        # b1, b2, b3 = principal_axes(bilayer)

        # Calculate angles between the various principal axes

        # The vector in the xy plane originating at the channel COM and extending to the COM of the TM2 vector
        # is normal to the plane against which we measure the TM2 tilting angle for radial tilt.
        for tm in tm2:
            a1, a2, a3 = principal_axes(tm)
            # Force the first axis to point up, always
            if a1[2] < 0: a1 = -a1
            radial_plane_normal = tm.centerOfMass() - channel.centerOfMass()
            radial_plane_normal[2] = 0
            radial_tilting_angle = 90-angle_in_degrees(a1, radial_plane_normal)

            # The plane relative to which we measure lateral tilting angle is orthogonal to the radial angle plane
            lateral_plane_normal = np.cross(radial_plane_normal, c1)
            lateral_tilting_angle = 90-angle_in_degrees(a1, lateral_plane_normal)

            print radial_tilting_angle, lateral_tilting_angle,

            # Print out diagnostic VMD arrow drawing commands if we see a really wacky lateral tilting angle
            if lateral_tilting_angle > 20:
                logging.info("Frame %d, radial=%.2f, lateral=%.2f" % (frame_num, radial_tilting_angle,
                                                                      lateral_tilting_angle))
                print "graphics top color blue"
                vmd_draw_arrow(channel.centerOfMass(), c1*20)
                vmd_draw_arrow(channel.centerOfMass(), radial_plane_normal)
                print "graphics top color red"
                vmd_draw_arrow(tm.centerOfMass(), a1*20)
                exit()

        print ""

    logging.info("Processed %d frames - all done." % frame_num)
