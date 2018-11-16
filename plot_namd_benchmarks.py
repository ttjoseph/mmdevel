#!/usr/bin/env python
#
# Uses matplotlib to plot a graph of ns/day as a function of number of cores, by parsing
# a bunch of NAMD logs
#
# Tom Joseph <thomas.joseph@uphs.upenn.edu>
import argparse
import re
import pprint
import numpy as np
import matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser(description='Use NAMD logs to make a performance plot')
    ap.add_argument('namdlogs', nargs='+', help='NAMD log file(s)')
    ap.add_argument('--title', default='NAMD scaling performance', help='Title of plot')
    ap.add_argument('--xlabel', default='Number of cores', help='Label for x axis')
    ap.add_argument('--ylabel', default='Nanoseconds per wallclock day', help='Label for y axis')
    ap.add_argument('--out', default='benchmark.png', help='Output filename')
    args = ap.parse_args()

    results = dict()

    for fname in args.namdlogs:
        lines = open(fname, 'r').readlines()
        # We want the most recent benchmark in each log
        for i in range(len(lines)-1, 0, -1):
            if 'days/ns' in lines[i]:
                num_cores = int(re.search(r'(\d+) CPUs', lines[i]).group(1))
                days_per_ns = float(re.search(r'([\.0-9]+) days/ns', lines[i]).group(1))

                if num_cores not in results:
                    results[num_cores] = []
                results[num_cores].append(1/days_per_ns)

                break

    pp = pprint.PrettyPrinter(indent=4)
    print('Nanoseconds per day by number of cores:')
    pp.pprint(results)

    avg, sd = [], []
    num_cores_list = list(sorted(results.keys()))
    for num_cores in num_cores_list:
        avg.append(np.mean(results[num_cores]))
        # We don't actually use standard deviation yet, but it's here if you
        # want to modify this script
        sd.append(np.std(results[num_cores]))

    plt.plot(num_cores_list, avg)
    plt.xticks(num_cores_list)
    plt.xlabel(args.xlabel or '')
    plt.ylabel(args.ylabel or '')
    plt.title(args.title or '')
    plt.savefig(args.out, dpi=300)



if __name__ == '__main__':
    main()