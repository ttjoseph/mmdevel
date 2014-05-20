#!/usr/bin/env python
#
# Shamelessly ripped off from:
#     http://www.ks.uiuc.edu/Research/vmd/script_library/scripts/orient/
import logging

from scipy import linalg
import MDAnalysis
import numpy as np

def inertia_tensor(atoms):
    """Somehow computes the inertia tensor, whatever that is. Does not mass weight."""
    # Remove center of mass
    coords = atoms.coordinates() - atoms.centerOfMass()
    ixx = ixy = ixz = iyy = iyz = izz = 0.0

    for c in coords:
        x, y, z = c[0], c[1], c[2]
#         set Ixx [expr $Ixx + $mm*($yy*$yy+$zz*$zz)]
        ixx += y*y+z*z
#         set Ixy [expr $Ixy - $mm*($xx*$yy)]
        ixy -= x*y
#         set Ixz [expr $Ixz - $mm*($xx*$zz)]
        ixz -= x*z
#         set Iyy [expr $Iyy + $mm*($xx*$xx+$zz*$zz)]
        iyy += x*x+z*z
#         set Iyz [expr $Iyz - $mm*($yy*$zz)]
        iyz -= y*z
#         set Izz [expr $Izz + $mm*($xx*$xx+$yy*$yy)]
        izz += x*x+y*y
#
#     return [list 2 3 3 $Ixx $Ixy $Ixz $Ixy $Iyy $Iyz $Ixz $Iyz $Izz]
    return np.array([[ixx, ixy, ixz], [ixy, iyy, iyz], [ixz, iyz, izz]])


def principal_axes(atoms):
    """Returns a matrix for which the columns are the principal axes."""
    it = inertia_tensor(atoms)
    # We can guarantee the inertia tensor is real and symmetric, so we can use eigh
    w, v = linalg.eigh(it)
    return v[:, 0], v[:, 1], v[:, 2]


def angle_in_degrees(u, v):
    """Returns the angle between two vectors, in degrees."""
    angle = np.arccos(np.dot(u, v)/np.linalg.norm(u)/np.linalg.norm(v))
    return angle*(180.0/np.pi)


if __name__ == '__main__':
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
    logging.info("Computes helix angles and whatever - Tom Joseph <thomas.joseph@mountsinai.org>")

    BasePath = '/Users/tom/Desktop/GLIC-CG/WithPFL/Intra-1/Restrained/'
    u = MDAnalysis.Universe(BasePath+'pr2.gro', BasePath+'md1-200kJ.xtc')
    # TM2 is resid 221:245
    # Atom indices below are valid only for the CG system.
    # TM2 A: bynum 510:562
    tm2 = u.selectAtoms('name BB and bynum 510:562')

    frame_num = 0
    for ts in u.trajectory:
        frame_num += 1
        a1, a2, a3 = principal_axes(tm2)
        z_axis = np.array([0, 0, 1])

        # As a first pass, calculate the angles with respect to the Z axis.
        # This represents the radial tilting angle.
        print frame_num, angle_in_degrees(a1, z_axis), angle_in_degrees(a2, z_axis), angle_in_degrees(a3, z_axis)

        # TODO: Lateral tilting angle. This is the angle against the plane formed by the Z axis and the vector
        # from the center of the pore toward the main principal axis of the helix
