#!/usr/bin/env python
#
# Generates a bunch of individual NAMD configurations representing lambda windows for
# a FEP calculation. Unlike embarassingly parallel FEP, this is only moderately
# embarassing, in that the transformation is broken into groups which are run serially.
# This way only the first lambda window of the group needs to be extensively equilibrated.
import numpy as np
import argparse
from os.path import isfile

sbatch_preambles = {}
sbatch_preambles['stampede'] = """#!/bin/bash
#SBATCH -J %(basename)s_%(group_id)d
#SBATCH -o %(basename)s_%(group_id)d.out
#SBATCH -e %(basename)s_%(group_id)d.err
#SBATCH -p normal         # queue
#SBATCH -N 4              # Number of nodes, not cores (16 cores/node)
#SBATCH -n 64             # Total number of MPI tasks (if omitted, n=N)
#SBATCH -t 48:00:00       # max time
#SBATCH --mail-user=thomas.joseph@uphs.upenn.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module restore system
module load gcc mvapich2

SRUN=ibrun
NAMD=~/bin/namd-2.12
"""

sbatch_preambles['stampede-knl'] = """#!/bin/bash
#SBATCH -J %(basename)s_%(group_id)d
#SBATCH -o %(basename)s_%(group_id)d.out
#SBATCH -e %(basename)s_%(group_id)d.err
#SBATCH -p normal         # queue
#SBATCH -N 2              # Number of nodes, not cores (68 real cores per node)
#SBATCH -n 8              # Total number of MPI tasks (if omitted, n=N)
#SBATCH -t 48:00:00       # max time
#SBATCH --mail-user=thomas.joseph@uphs.upenn.edu
#SBATCH --mail-type=begin --mail-type=end

module restore system
module load intel/16.0.3
module load impi/5.1.3
module load namd/2.12

SRUN=ibrun
NAMD="$TACC_NAMD_BIN/namd2 +ppn 32 +pemap 0-63+68 +commap 64-67"
"""

sbatch_preambles['perceval'] = """#!/bin/bash
#SBATCH -J %(basename)s_%(group_id)d
#SBATCH -o %(basename)s_%(group_id)d.out
#SBATCH -e %(basename)s_%(group_id)d.err
#SBATCH -p main
#SBATCH -x memnode001 # Don't use this node as it is real slow
#SBATCH -N 1 -n 24
#SBATCH -t 48:00:00
#SBATCH --mem=12000
#SBATCH --export=ALL

# Run these commands before submitting this script:
#   module purge
#   module load gcc mvapich2

SRUN="srun -n $SLURM_NTASKS --mpi=pmi2"
NAMD=$HOME/bin/namd2-2.12

echo "Starting %(basename)s_%(group_id)d" | ~/bin/slack-ttjoseph

"""



ap = argparse.ArgumentParser(description='Generate mildly embarrasingly parallel FEP input for NAMD')
ap.add_argument('basename', help='Base filename, e.g. "fep", such that we generate fep<n>.namd which sources common_fep.namd')
ap.add_argument('inputname', help='Per-group input system name (e.g. equilibrated system)')
ap.add_argument('clustername', help='Supercomputer cluster you are using', default='perceval')
ap.add_argument('--num-windows', type=int, default=30, help='Number of lambda windows')
ap.add_argument('--windows-per-group', type=int, default=5, help='Maximum number of lambda windows per batch job')
ap.add_argument('--per-group-num-equil-steps', type=int, default=500000)
ap.add_argument('--per-window-num-equil-steps', type=int, default=200000)
ap.add_argument('--per-window-prod-steps', type=int, default=1000000)
ap.add_argument('--continue', dest='continue_prod', action='store_true', help='Continue from previous run as specified by inputname')
ap.set_defaults(continue_prod=False)
args = ap.parse_args()

x = np.linspace(0, 1.0, num=(args.num_windows+1), endpoint=True)
all_lambdas = np.sin((np.pi/2)*x)
print('Here are the %d lambda windows:' % args.num_windows)
print(all_lambdas)
group_id = 0

for i in range(len(all_lambdas)-1):
    starting_a_group = True if i % args.windows_per_group == 0 else False
    if starting_a_group:
        group_id += 1

    data = {
        'basename': args.basename,
        'inputname': args.inputname,
        'outputname': '%s%03d' % (args.basename, i),
        'group_id': group_id,
        'l0': all_lambdas[i],
        'l1': all_lambdas[i+1]
    }
    if starting_a_group:
        data['alchequilsteps'] = args.per_group_num_equil_steps
    else:
        data['inputname'] = '%s%03d' % (args.basename, i - 1)
        data['alchequilsteps'] = args.per_window_num_equil_steps

    # If we are adding on more production sampling to a previous run, we don't need to
    # do any equilibration, and we should start from the same window from the previous run
    if args.continue_prod is True:
        data['alchequilsteps'] = 0
        data['inputname'] = '%s%03d' % (args.inputname, i)

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
        s = sbatch_preambles[args.clustername] % data
        with open('submit-%s-group-%d.sh' % (args.basename, group_id), 'w') as f:
            f.write(s)
            for j in range(i, i+args.windows_per_group):
                f.write('$SRUN $NAMD %s%03d.namd > %s%03d.log\n' % (args.basename, j, args.basename, j))
                f.write('echo "%s%03d (group %d) done" | ~/bin/slack-ttjoseph\n\n' % (args.basename, j, group_id))
