#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import sys

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Calculate descriptive statistics for columns of numbers in a CSV')
    ap.add_argument('csvfile', help='CSV filename')
    args = ap.parse_args()

    d = pd.read_csv(args.csvfile)
    print 'Label,Mean,SD'
    for col in d:
        print '%s,%f,%f' % (col, np.mean(d[col]), np.std(d[col]))
