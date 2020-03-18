#!/usr/bin/env python
#
# Generate a plot of a CSV file with a bunch of columns.
# Especially useful for time-series data, such as RMSD over an MD trajectory
import argparse
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


# Draw broken bars as an annotation. This was intended to help draw protein domains.
def do_broken_barh(csvfile, plt):
    df = pd.read_csv(csvfile)
    # Each horizontal bar is part of a group with a name.
    groups = sorted(list(df['group'].unique()))
    y_bottom, y_top = plt.gca().get_ylim()
    counter = 0
    colors = ['blue', 'red', 'green', 'purple', 'orange']

    print("do_broken_barh: All sorts of hackiness in here. You might want to edit to source to get the output you want.")
    y_top -= 3
    for group in groups:
        rows = df.query(f'group == "{group}"').sort_values(by=['xmin'])
        print(rows)
        plt.broken_barh(list(zip(rows['xmin'], rows['xwidth'])), (y_top-counter, 0.8), facecolors=colors[counter])
        for idx, row in rows.iterrows():
            plt.text(row['xmin'] + (row['xwidth'] / 2), y_top-counter+0.4, row['title'],
                horizontalalignment='center', verticalalignment='center', fontsize='x-small',
                color='white')
        counter += 1


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Plot a CSV of time-series data')
    ap.add_argument('in_csv', help='Input CSV file')
    ap.add_argument('out_img', help='Output image file')
    ap.add_argument('--columns', help='Column indices to include, comma-separated, starting from 0 (default: all)')
    ap.add_argument('--x-label', default='', help='X-axis label')
    ap.add_argument('--y-label', default='', help='Y-axis label')
    ap.add_argument('--x-scale', type=float, default=1, help='Scale index of X axis data by this much (e.g. to convert trajectory frame number to time)')
    ap.add_argument('--index-offsets', help='Add these numbers, comma-separated, to each respective index value, before scaling with --x-scale')
    ap.add_argument('--legend-labels', help='Comma-separated labels for the legend. Yes, that means you can\'t include a comma in a label')
    ap.add_argument('--broken-barh-csv', help='CSV file specifying broken horizontal bars. Coordinates are post-transformation. Look at my source code')
    args = ap.parse_args()

    df = pd.read_csv(args.in_csv)

    # Pull out only the columns the user specifies, if applicable
    if args.columns is not None:
        cols = [int(x) for x in args.columns.split(',')]
        df = df.iloc[:, cols]

    # Extract user-specified scale and offset values to apply to index
    index_offsets = [int(x) for x in args.index_offsets.split(',')] if args.index_offsets is not None else np.zeros(len(df.columns))
    index_scale = np.array(df.index)*np.array(args.x_scale)

    if len(index_offsets) != len(df.columns):
        print(f"You specified {len(index_offsets)} index offsets but you have {len(df.columns)} columns in your data.")

    mpl.rcParams['lines.linewidth'] = 0.25

    print(df.columns)

    i = 0
    for col_name, col_data in df.iteritems():
        this_index = (np.array(df.index) + index_offsets[i])*np.array(args.x_scale)
        plt.plot(this_index, col_data)
        i += 1

    plt.xlabel(args.x_label)
    plt.ylabel(args.y_label)
    labels = None
    if args.legend_labels is not None:
        labels = args.legend_labels.split(',')

    # Do broken bar stuff
    if args.broken_barh_csv is not None:
        do_broken_barh(args.broken_barh_csv, plt)

    plt.legend(fontsize='x-small', loc='upper right', labels=labels)
    if labels == ['no']:
        plt.gca().get_legend().remove()
    plt.gcf().set_size_inches(7,4)
    plt.savefig(args.out_img, dpi=600)
