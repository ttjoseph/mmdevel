#!/usr/bin/env python
#
# Generate a plot of a CSV file with a bunch of columns.
# Especially useful for time-series data, such as RMSD over an MD trajectory
import argparse
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Plot a CSV of time-series data')
    ap.add_argument('in_csv', help='Input CSV file')
    ap.add_argument('out_img', help='Output image file')
    ap.add_argument('--columns', help='Column indices to include, comma-separated, starting from 0 (default: all)')
    ap.add_argument('--x-label', default='', help='X-axis label')
    ap.add_argument('--y-label', default='', help='Y-axis label')
    ap.add_argument('--x-scale', type=float, default=1, help='Scale index of X axis data by this much (e.g. to convert trajectory frame number to time)')
    ap.add_argument('--legend-labels', help='Comma-separated labels for the legend. Yes, that means you can\'t include a comma in a label')
    args = ap.parse_args()

    df = pd.read_csv(args.in_csv)

    if args.columns is not None:
        cols = [int(x) for x in args.columns.split(',')]
        df = df.iloc[:, cols]

    df.set_index(np.array(df.index)*np.array(args.x_scale))

    mpl.rcParams['lines.linewidth'] = 0.25

    print(df.columns)

    df.plot()
    plt.xlabel(args.x_label)
    plt.ylabel(args.y_label)
    labels = None
    if args.legend_labels is not None:
        labels = args.legend_labels.split(',')

    plt.legend(fontsize='x-small', loc='upper right', labels=labels)
    if labels == ['no']:
        plt.gca().get_legend().remove()
    plt.gcf().set_size_inches(7,4)
    plt.savefig(args.out_img, dpi=600)
