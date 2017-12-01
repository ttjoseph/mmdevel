import argparse
import re
import sys
import os
import subprocess
import numpy as np
import MDAnalysis as mda
# From the VMD documentation of the DX plugin
#   Format of the file is:
#    # Comments
#    .
#    .
#    .
#    object 1 class gridpositions counts xn yn zn
#    origin xorg yorg zorg
#    delta xdel 0 0
#    delta 0 ydel 0
#    delta 0 0 zdel
#    object 2 class gridconnections counts xn yn zn
#    object 3 class array type double rank 0 items [ xn*yn*zn ] data follows
#    f1 f2 f3
#    f4 f5 f6
#    .
#    .
#    .
#    object "Dataset name" class field
#
#    Where xn, yn, and zn are the number of data points along each axis;
#    xorg, yorg, and zorg is the origin of the grid, in angstroms;
#    xdel, ydel, and zdel are the scaling factors to convert grid units to
#    angstroms.
#    Grid data follows, with three values per line, ordered z fast, y medium,
#    and x slow.

class ElectrostaticPotential(object):
    """Encapsulates a .dx file which probably contains a electrostatic potential field.

    The use case (for now) is to extract the mean electrostatic potential for a subset of
    points. For example, the mean electrostatic potential in the solvent only. For this
    the caller could pass in a list of coordinates (in Angstroms) and we could do the
    averaging or whatever and return a scalar in units of say kT/e.

    Only meant to deal with APBS style DX files, and only meant to be as good as VMD's
    parser, which is a pretty low standard."""

    def __init__(self):
        self.rawgrid = None

    def __init__(self, dx_fname):
        """Load the thing from a .dx file"""
        f = open(dx_fname)
        # Parse comments
        line = '#'
        while line.startswith('#'):
            line = f.readline()

        # object 1 class gridpositions counts xn yn zn
        m = re.match('object 1 class gridpositions counts (\d+) (\d+) (\d+)', line)
        self.num_points = int(m.group(1)), int(m.group(2)), int(m.group(3))
        # origin xorg yorg zorg
        m = re.match('origin ([-\d\.]+) ([-\d\.]+) ([-\d\.]+)', f.readline())
        self.origin = float(m.group(1)), float(m.group(2)), float(m.group(3))
        # Get delta values
        # delta xdel 0 0
        # delta 0 ydel 0
        # delta 0 0 zdel
        self.delta = (float(re.match('delta ([-\d\.]+) [-\d\.]+ [-\d\.]+', f.readline()).group(1)),
                      float(re.match('delta [-\d\.]+ ([-\d\.]+) [-\d\.]+', f.readline()).group(1)),
                      float(re.match('delta [-\d\.]+ [-\d\.]+ ([-\d\.]+)', f.readline()).group(1)))

        # object 2 class gridconnections counts xn yn zn
        m = re.match('object 2 class gridconnections counts (\d+) (\d+) (\d+)', f.readline())
        self.num_points = np.array([int(m.group(1)), int(m.group(2)), int(m.group(3))])
        # object 3 class array type double rank 0 items [ xn*yn*zn ] data follows
        num_vals = int(re.match('object 3 class array type double rank 0 items (\d+) data follows',
                                f.readline()).group(1))

        # print 'Origin:', self.origin
        # print 'Counts:', self.num_points
        # print 'Deltas:', self.delta
        assert num_vals == np.prod(self.num_points)

        # Suck in the actual data
        data = np.zeros(num_vals)
        counter = 0
        while True:
            line = f.readline()
            if line.startswith('attr') or line.startswith('obj'):
                break
            vals = [float(x) for x in line.split()]
            for x in vals:
                data[counter] = x
                counter += 1

        if counter != num_vals:
            print >>sys.stderr, 'Uh oh. I read %d values but there were supposed to be %d.' % (counter, num_vals)

        self.rawgrid = data.reshape(self.num_points)

    # This is useful when for whatever reason the entire potential is translated relative to
    # the effective origin of the protein coordinates. Still don't know why this happens.
    def set_origin(self, origin):
        self.origin = origin

    def mean_potential_serial(self, points):
        vals = np.zeros((len(points), 3))
        counter = 0
        for p in points:
            # Get nearest coordinate we have
            grid_point = list(int(x) for x in (p - self.origin) / self.delta)
            vals[counter] = self.rawgrid[grid_point[0], grid_point[1], grid_point[2]]
            counter += 1
           
        return np.mean(vals)

    def mean_potential(self, points):
        def nearest_val(p):
            gp = np.clip((p - self.origin) / self.delta, (0, 0, 0), self.num_points - 1)
            return self.rawgrid[tuple(gp.astype(int))]

        return np.mean(np.apply_along_axis(nearest_val, 1, points))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Calculate electrostatic potential correction necessary for FEP/PME/charged ligand')
    ap.add_argument('psf', help='NAMD PSF topology file')
    ap.add_argument('dcd', help='Trajectory file')
    ap.add_argument('--vmdbin', help='VMD binary', default='vmd')
    ap.add_argument('--keep-temp-files', help='Keep temporary files such as frames_*.dx and frames.dcd', action='store_true')
    ap.set_defaults(keep_temp_files=False)
    args = ap.parse_args()

    print >>sys.stderr, 'Writing PBC-wrapped trajectory to frames.dcd and computing PME potentials...'
    p = subprocess.Popen([args.vmdbin, '-dispdev', 'text'], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate("""package require pmepot
    package require pbctools
    mol new %s
    mol addfile %s waitfor all
    set numframes [molinfo top get numframes]
    pbc wrap -all
    animate write dcd frames.dcd
    for {set i 0} {$i < $numframes} {incr i} {
      pmepot -frames $i:$i -loadmol none -dxfile [format "frame_%%d.dx" $i]
    }
    quit""" % (args.psf, args.dcd))

    print >>sys.stderr, 'Loading %s and %s.' % (args.psf, 'frames.dcd')
    u = mda.Universe(args.psf, 'frames.dcd')
    # The atom selection language doesn't seem to have a VMD-like "within" keyword so we have to do this bullshit
    print >>sys.stderr, 'Selecting the atoms.'
    solvent = u.select_atoms('not (around 5 (protein or resname POPC CHL1)) and not (protein or resname POPC CHL1)',
                             updating=True)
    print >>sys.stderr, 'Selected %d atoms.' % len(solvent.atoms)

    counter, potentials = 0, np.zeros(u.trajectory.n_frames)
    for ts in u.trajectory:
        ep = ElectrostaticPotential('frame_%d.dx' % counter)
        minx, miny, minz = np.min(u.atoms.positions, axis=0)
        print >>sys.stderr, 'Minimum coordinates in DCD file, which we ignore:', (minx, miny, minz)
        print >>sys.stderr, 'Origin in DX file:', ep.origin
        potentials[counter] = ep.mean_potential(solvent.positions)
        print 'Frame %d mean solvent potential is: %.2f kT/e' % (counter, potentials[counter])
        if args.keep_temp_files is False:
            os.unlink('frame_%d.dx' % counter)
        counter += 1

    print 'Mean solvent potential across all %d frames: %.2f kT/e' % (u.trajectory.n_frames, np.mean(potentials))
    if args.keep_temp_files is False:
        os.unlink('frames.dcd')