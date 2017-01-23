#!/usr/bin/env python
#
# Generates a bunch of individual NAMD configurations representing lambda windows for
# a FEP calculation. Unlike embarassingly parallel FEP, this is only moderately
# embarassing, in that the transformation is broken into groups which are run serially.
# This way only the first lambda window of the group needs to be extensively equilibrated.
import numpy as np
import argparse
from os.path import isfile

ap = argparse.ArgumentParser(description='Generate mildly embarrasingly parallel FEP input for NAMD')
ap.add_argument('basename', help='Base filename, e.g. "fep", such that we generate fep<n>.namd which sources common_fep.namd')
ap.add_argument('inputname', help='Per-group input system name (e.g. equilibrated system)')
ap.add_argument('--more-steps', action='store_true', help='Actually, do more sampling for existing FEP calculation, and interpret inputname as the (numeric) label')
ap.add_argument('--num-windows', default=50, help='Number of lambda windows')
ap.add_argument('--num-groups', default=10, help='Number of groups. Must be evenly divisible into num-windows.')
ap.add_argument('--per-group-num-equil-steps', default=500000)
ap.add_argument('--per-window-num-equil-steps', default=200000)
ap.add_argument('--per-window-prod-steps', default=1000000)
args = ap.parse_args()

if args.num_windows % args.num_groups != 0:
    print('Cannot evenly divide %d windows into %d groups.' % (args.num_windows, args.num_groups))

if args.more_steps:
    print('You are adding on more production sampling. So I will assume you do not want more equilibration.')
    # Ensure inputname can be cleanly converted to an integer
    # If it's 0, then assume we are starting from the unlabeled first run
    try:
        int(args.inputname)
    except:
        print('But in this case, you need to provide an integer for the inputname, which will be used as a suffix')
        print('to the basename. If this is the first time you are doing --more-steps, just say 0.')
        exit(1)

    # Ensure that .restart.coor, .restart.vel, .restart.xsc files that we'll need actually exist.
    # This is to protect against user typos.
    for i in range(args.num_windows):
        for suffix in ('coor', 'vel', 'xsc'):
            if int(args.inputname) > 0:
                fname = '%s%03d_%d.restart.%s' % (args.basename, i, int(args.inputname), suffix)
            else:
                fname = '%s%03d.restart.%s' % (args.basename, i, suffix)
            if isfile(fname) is False:
                print('%s does not exist, and you need it. Probably other stuff is wrong too!' % fname)
                exit(1)

group_size = args.num_windows / args.num_groups
all_lambdas = np.linspace(0, 1.0, num=50, endpoint=False)
delta = all_lambdas[1] - all_lambdas[0]
for i in range(len(all_lambdas)):
    starting_a_group = True if i % group_size == 0 else False
    group_id = i / group_size

    l0, l1 = all_lambdas[i], all_lambdas[i] + delta
    data = {
        'basename': args.basename,
        'inputname': args.inputname,
        'outputname': '%s%03d' % (args.basename, i),
        'group_id': group_id,
        'l0': l0,
        'l1': l1
    }
    if starting_a_group:
        data['alchequilsteps'] = args.per_group_num_equil_steps
    else:
        data['inputname'] = '%s%03d' % (args.basename, i - 1)
        data['alchequilsteps'] = args.per_window_num_equil_steps

    if args.more_steps:
        if int(args.inputname) == 0:
            data['inputname'] = '%s%03d.restart' % (args.basename, i)
        else:
            data['inputname'] = '%s%03d_%d.restart' % (args.basename, i, int(args.inputname))

        data['outputname'] = '%s%03d_%d' % (args.basename, i, int(args.inputname) + 1)
        data['alchequilsteps'] = 0

    data['totalsteps'] = data['alchequilsteps'] + args.per_window_prod_steps
    s = """# Do FEP for lambda %(l0)f to %(l1)f
set inputname %(inputname)s;
set outputname %(outputname)s;
alchLambda %(l0)f
alchLambda2 %(l1)f
source common_%(basename)s.namd
alchEquilSteps %(alchequilsteps)d
run %(totalsteps)d
""" % data
    with open('%s.namd' % data['outputname'], 'w') as f:
        f.write(s) 

    # Make a submit script for this window group, for convenience
    if starting_a_group:
        s = """#!/bin/bash
#SBATCH -J %(basename)s_%(group_id)d
#SBATCH -o %(basename)s_%(group_id)d.out
#SBATCH -e %(basename)s_%(group_id)d.err
#SBATCH -p normal         # queue
#SBATCH -N 10             # Number of nodes, not cores (16 cores/node)
#SBATCH -n 160            # Total number of MPI tasks (if omitted, n=N)
#SBATCH -t 48:00:00       # max time
#SBATCH --mail-user=thomas.joseph@uphs.upenn.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module restore system
module load gcc mvapich2

""" % data
        with open('submit-%s-group-%d.sh' % (args.basename, group_id), 'w') as f:
            f.write(s)
            for j in range(i, i+group_size):
                suffix = '_%d' % (int(args.inputname) + 1) if args.more_steps else ''
                f.write('ibrun ~/bin/namd2-2.12 %s%03d%s.namd > %s%03d%s.log\n' % (args.basename, j, suffix, args.basename, j, suffix))
                f.write('echo "%s%03d%s (group %d) done" | ~/bin/slack-ttjoseph\n\n' % (args.basename, j, suffix, group_id))
