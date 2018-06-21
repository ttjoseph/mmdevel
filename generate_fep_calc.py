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
ap.add_argument('basenamdconf', help='Base filename, e.g. "fep", such that we generate fep<n>.namd which sources common_fep.namd')
ap.add_argument('startsystem', help='Per-group input system name (e.g. equilibrated system)')
ap.add_argument('--num-windows', type=int, default=30, help='Number of lambda windows')
ap.add_argument('--left-dense', dest='left_dense', action='store_true', help='Denser lambda coverage closer to 0, rather than 1?')
ap.add_argument('--idws', dest='idws', action='store_true', help='Use interleaved double-wide sampling')
ap.add_argument('--windows-per-group', type=int, default=5, help='Maximum number of lambda windows per batch job')
ap.add_argument('--per-group-num-equil-steps', type=int, default=500000)
ap.add_argument('--per-window-num-equil-steps', type=int, default=200000)
ap.add_argument('--per-window-prod-steps', type=int, default=1000000)
ap.add_argument('--continue', dest='continue_prod', action='store_true', help='Continue from previous run as specified by startsystem')
ap.set_defaults(continue_prod=False, left_dense=False, idws=False)
args = ap.parse_args()

x = np.linspace(0, 1.0, num=(args.num_windows+1), endpoint=True)
# Make the lambdas denser to one side or the other of the FEP run.
# Often, there will be big dG changes with low or high lambda values, so this is intended
# to sample those areas more by using smaller windows.
all_lambdas = np.sin((np.pi/2)*x) if args.left_dense is False else (1-np.cos((np.pi/2)*x))
print('Here are the %d lambda windows:' % args.num_windows)
print(all_lambdas)
group_id = 0

for i in range(len(all_lambdas)-1):
    starting_a_group = True if i % args.windows_per_group == 0 else False
    if starting_a_group:
        group_id += 1

    data = {
        'basenamdconf': args.basenamdconf,
        'startsystem': args.startsystem,
        'outputname': '%s%03d' % (args.basenamdconf, i),
        'group_id': group_id,
        'l0': all_lambdas[i],
        'l1': all_lambdas[i+1],
        'idws_string': 'alchLambdaIDWS %f' % all_lambdas[i-1] if i > 0 else ''
    }
    if starting_a_group:
        data['alchequilsteps'] = args.per_group_num_equil_steps
    else:
        data['startsystem'] = '%s%03d' % (args.basenamdconf, i - 1)
        data['alchequilsteps'] = args.per_window_num_equil_steps

    # If we are adding on more production sampling to a previous run, we don't need to
    # do any equilibration, and we should start from the same window from the previous run
    if args.continue_prod is True:
        data['alchequilsteps'] = 0
        data['startsystem'] = '%s%03d' % (args.startsystem, i)

    data['totalsteps'] = data['alchequilsteps'] + args.per_window_prod_steps
    s = """# Do FEP for lambda %(l0)f to %(l1)f
set inputname %(startsystem)s;
set outputname %(outputname)s;
alchLambda %(l0)f
alchLambda2 %(l1)f
%(idws_string)s
source common_%(basenamdconf)s.namd
alchEquilSteps %(alchequilsteps)d
run %(totalsteps)d
""" % data
    with open('%s.namd' % data['outputname'], 'w') as f:
        f.write(s)
