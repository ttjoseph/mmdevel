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
args = ap.parse_args()

# Get latest prod?.restart.coor, sorted properly numerically
# Because we don't want to pull in natsort or similar, lex-sort by zero padded number
skel_filename = '%s.skel.namd' % args.prefix
if os.path.isfile(skel_filename) is False:
    exit('Skeleton NAMD config file %s does not seem to be there' % skel_filename)
files = sorted(glob.glob('%s*.restart.coor' % args.prefix), key=lambda k: '%09d' % int(re.search('\d+', k).group()))
if len(files) == 0:
    exit('Unable to find any NAMD restart files')
last_restart_filename = files[-1]
# Extract the last index number from the last restart filename
last_index = int(re.search('\d+', last_restart_filename).group())
# Next index should be one higher than the last one
next_index = last_index + 1
next_config_filename = '%s%d.namd' % (args.prefix, next_index)
# Create NAMD config file from the last config file and setting variables
with open(next_config_filename, 'w') as next_config_f:
    next_config_f.write("""set inputname %s%d.restart;
set outputname %s%d;
### Skel file below here ###
""" % (args.prefix, last_index, args.prefix, next_index))
    with open(skel_filename) as skel_f:
        next_config_f.write(skel_f.read())

# Caller shell script is intended to capture this output
print(next_config_filename)