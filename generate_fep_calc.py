#!/usr/bin/env python
#
# Generates a bunch of individual NAMD configurations representing lambda windows for
# a FEP calculation. Unlike embarassingly parallel FEP, this is only moderately
# embarassing, in that the transformation is broken into groups which are run serially.
# This way only the first lambda window of the group needs to be extensively equilibrated.
import numpy as np
from scipy.special import erf
import argparse
from os.path import isfile

ap = argparse.ArgumentParser(description='Generate mildly embarrasingly parallel FEP input for NAMD')
ap.add_argument('basenamdconf', help='Base filename, e.g. "fep", such that we generate fep<n>.namd which sources common_fep.namd')
ap.add_argument('startsystem', help='Per-group input system name (e.g. equilibrated system)')
ap.add_argument('--num-windows', type=int, default=30, help='Number of lambda windows')
# ap.add_argument('--left-dense', dest='left_dense', action='store_true', help='Denser lambda coverage closer to 0, rather than 1?')
ap.add_argument('--lambda-scheme', '-s', default='erf-5', help='Lambda schedule scheme')
ap.add_argument('--no-idws', dest='idws', action='store_false', help='Do not use interleaved double-wide sampling')
ap.add_argument('--minimize', dest='minimize', action='store_true', help='Minimize at the start of each group')
ap.add_argument('--windows-per-group', type=int, default=5, help='Maximum number of lambda windows per batch job')
ap.add_argument('--per-group-num-equil-steps', type=int, default=500000)
ap.add_argument('--per-window-num-equil-steps', type=int, default=200000)
ap.add_argument('--per-window-prod-steps', type=int, default=1000000)
ap.add_argument('--continue', dest='continue_prod', action='store_true', help='Continue from previous run as specified by startsystem')
ap.set_defaults(continue_prod=False, left_dense=False, idws=True)
args = ap.parse_args()

x = np.linspace(0, 1.0, num=(args.num_windows+1), endpoint=True)
# Make the lambdas denser to one side or the other of the FEP run.
# Often, there will be big dG changes with low or high lambda values, so this is intended
# to sample those areas more by using smaller windows.
if args.lambda_scheme.startswith('erf-'):
    # \frac{1}{2}\operatorname{erf}\left(4\left(x-\frac{1}{\left(2\right)}\right)\right)+\frac{1}{2}
    coeff = float(args.lambda_scheme[4:])
    all_lambdas = 0.5 * erf(coeff * (x - 0.5)) + 0.5
    all_lambdas[0], all_lambdas[-1] = 0.0, 1.0
elif args.lambda_scheme == 'right-dense-cosine':
    all_lambdas = np.sin((np.pi/2)*x)
elif args.lambda_scheme == 'left-dense-cosine':
    all_lambdas = 1-np.cos((np.pi/2)*x)
elif args.lambda_scheme == 'right-dense-log':
    all_lambdas = 1-np.log(np.power(x, -1/3))
elif args.lambda_scheme == 'left-dense-log':
    all_lambdas = 1-np.log(np.power(1-x, -1/3))
elif args.lambda_scheme == 'evenly-spaced':
    all_lambdas = x
else:
    print(f"Error: Somehow, you specified lambda scheme {args.lambda_scheme}, which I don't know how to do.")
    exit(1)
print('Here are the %d lambdas:' % args.num_windows)
np.set_printoptions(precision=5, suppress=True)
print(all_lambdas)
real_num_windows = args.num_windows if args.idws is False else args.num_windows+1
is_last_idws_window = False

print('That means %d windows.' % real_num_windows)

for i in range(real_num_windows):
    starting_a_group = True if i % args.windows_per_group == 0 else False

    is_last_idws_window = (args.idws is True and i == (real_num_windows-1))

    l0, l1 = None, None

    if is_last_idws_window:
        l0, l1 = all_lambdas[i], all_lambdas[i-1]
    else:
        l0, l1 = all_lambdas[i], all_lambdas[i+1]

    data = {
        'basenamdconf': args.basenamdconf,
        'startsystem': args.startsystem,
        'outputname': '%s%03d' % (args.basenamdconf, i),
        'l0': l0,
        'l1': l1,
        'idws_string': 'alchLambdaIDWS %f' % all_lambdas[i-1] if i > 0 and is_last_idws_window is False else '',
        'minimize_string': ''
    }
    if starting_a_group:
        data['alchequilsteps'] = args.per_group_num_equil_steps
        # Minimize at the start of every group
        data['minimize_string'] = 'minimize 1000' if args.minimize is True else ''
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
# We must set first time step 0 right before running dynamics
# If it isn't, the fepout format will be screwed up
%(minimize_string)s
firsttimestep 0
alchEquilSteps %(alchequilsteps)d
run %(totalsteps)d
""" % data
    with open('%s.namd' % data['outputname'], 'w') as f:
        f.write(s)
