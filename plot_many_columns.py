#!/usr/bin/env python
#
# Generate a plot of a CSV file with a bunch of columns.
# Especially useful for time-series data, such as RMSD over an MD trajectory
import argparse
import pandas as pd
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

    mpl.rcParams['lines.linewidth'] = 0.25
    
    df.plot(x=df.index*args.x_scale, y=df.columns)
    plt.xlabel(unicode(args.x_label, 'utf8'))
    plt.ylabel(unicode(args.y_label, 'utf8'))
    labels = None
    if args.legend_labels is not None:
        labels = args.legend_labels.split(',')

    plt.legend(fontsize='x-small', loc='upper right', labels=labels)
    if labels == ['no']:
        plt.gca().get_legend().remove()
    plt.gcf().set_size_inches(7,4)
    plt.savefig(args.out_img, dpi=600)
