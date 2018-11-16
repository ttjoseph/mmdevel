#!/usr/bin/env python
#
# Uses matplotlib to make a histogram from data in a CSV file.
import argparse
import sys
from os.path import splitext, basename
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Make a histogram figure out of CSV file columns')
    ap.add_argument('csvfile', help='CSV file that contains a header with labels')
    ap.add_argument('--title', help='Title of plot')
    ap.add_argument('--xlabel', help='Label for x axis')
    ap.add_argument('--ylabel', help='Label for y axis')
    args = ap.parse_args()

    df = pd.read_csv(args.csvfile, header=0)

    plt.figure()

    for column in df:
        plt.hist(df[column], bins=50, label=column)

    plt.xlabel(args.xlabel or '')
    plt.ylabel(args.ylabel or '')
    plt.title(args.title or '')
    plt.legend()
    out_filename = '%s-hist.png' % splitext(basename(args.csvfile))[0]
    plt.savefig(out_filename)
    print >>sys.stderr, 'Wrote histogram to %s.' % out_filename
