#!/usr/bin/python
#
# Creates config file for "next" NAMD simulation.
# For example, if you have prod2.restart.coor, then it will use prod.skel.namd to
# automatically create a prod3.namd
import glob
import re
import os
from sys import exit
import argparse

ap = argparse.ArgumentParser(description='Create and print the "next" NAMD production simulation')
ap.add_argument('--prefix', default='prod', help='Prefix to NAMD files: e.g. "prod"')
ap.add_argument('--inputname', '-i', help='Input .coor file to use if no "prod" restart files available')
args = ap.parse_args()

# Get latest prod?.restart.coor, sorted properly numerically
# Because we don't want to pull in natsort or similar, lex-sort by zero padded number
skel_filename = '%s.skel.namd' % args.prefix
if os.path.isfile(skel_filename) is False:
    exit('Skeleton NAMD config file %s does not seem to be there' % skel_filename)
files = sorted(glob.glob('%s*.restart.coor' % args.prefix), key=lambda k: '%09d' % int(re.search(r'\d+', k).group()))

# No production restart file? Start from the end of equilibration
start_from_equilibration = False
# Maybe the user did some other min/eq protocol, or renamed the CHARMM-GUI ones
# If they specified something explicit for inputname we should fail if we can't find it,
# rather than falling back to what is probably the wrong one.
if args.inputname is None:
    inputnames_to_glob = ['step6.6_equilibration.restart.coor', 'step4_equilibration.restart.coor', 'mineq.restart.coor']
else:
    inputnames_to_glob = [args.inputname,]

if len(files) == 0:
    for prev_run in inputnames_to_glob:
        files = glob.glob(prev_run)
        if len(files) > 0:
            start_from_equilibration = True
            break

if len(files) == 0:
    exit(f'Unable to find any NAMD restart files. I looked for {args.prefix}*.restart.coor and {inputnames_to_glob}')
last_restart_filename = files[-1]
last_index = 0

# Extract the last index number from the last restart filename
# If we are starting from equilibration, suppress fancy prediction of filenames
if start_from_equilibration:
    inputname = last_restart_filename.replace('.restart.coor', '')
else:
    last_index = int(re.search('\d+', last_restart_filename).group())
    inputname = '%s%d.restart' % (args.prefix, last_index)

# Next index should be one higher than the last one
next_index = last_index + 1
next_config_prefix = '%s%d' % (args.prefix, next_index)

# Create NAMD config file from the last config file and setting variables
with open('%s.namd' % next_config_prefix, 'w') as next_config_f:
    next_config_f.write("""set inputname %s;
set outputname %s%d;
### Skel file below here ###
""" % (inputname, args.prefix, next_index))
    with open(skel_filename) as skel_f:
        next_config_f.write(skel_f.read())

# Caller shell script is intended to capture this output
print(next_config_prefix)
