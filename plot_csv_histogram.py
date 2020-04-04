#!/usr/bin/env python3
#
# Uses matplotlib to make a histogram from data in a CSV file.
import argparse
import sys
from os.path import splitext, basename
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Plot time series from CSV, with time in first column')
    ap.add_argument('csvfile', help='CSV file that contains a header with labels')
    ap.add_argument('--title', help='Title of plot')
    ap.add_argument('--xlabel', help='Label for x axis')
    ap.add_argument('--ylabel', help='Label for y axis')
    ap.add_argument('--columns', help='Keep only these column numbers')
    args = ap.parse_args()

    df = pd.read_csv(args.csvfile, header=0)

    plt.figure()

    for column in df:
        plt.hist(df[column], bins=50, label=column)

    plt.xlabel(args.xlabel or '')
    plt.ylabel(args.ylabel or '')
    plt.title(args.title or '')
    plt.legend()
    out_filename = '{}-hist.png'.format(splitext(basename(args.csvfile))[0])
    plt.savefig(out_filename)
    print('Wrote histogram to {}.'.format(out_filename), file=sys.stderr)
