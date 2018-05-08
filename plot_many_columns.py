#!/usr/bin/env python
import argparse
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Plot a CSV of time-series data')
    ap.add_argument('in_csv', help='Input CSV file')
    ap.add_argument('out_pdf', help='Output PDF file')
    args = ap.parse_args()

    df = pd.read_csv(args.in_csv)
    mpl.rcParams['lines.linewidth'] = 1
    plt.figure(figsize=(7.5, 5))
    df.plot(x=df.index, y=df.columns)
    plt.savefig(args.out_pdf)
