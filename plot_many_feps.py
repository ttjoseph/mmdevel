#!/usr/bin/env python3
import argparse
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_fepout(fname):
	with open(fname) as f:
		lines = f.readlines()
		fepenergy, deltas, lambdas = list(), list(), list()
		for line in lines:
			tokens = line.split()
			if line.startswith('FepEnergy'):
				# A single sample
				# The energy is always in the tenth column
				fepenergy.append(float(tokens[9]))
			elif line.startswith('#Free'):
				# delta-G at end of window, in the twelfth column
				deltas.append(float(tokens[11]))
			elif line.startswith('#NEW FEP WINDOW'):
				# Lambda value in seventh column
				lambdas.append(float(tokens[6]))

		return fepenergy, deltas, lambdas


def main():
	ap = argparse.ArgumentParser(description='Make a big plot of FEP curves')
	ap.add_argument('dirname', nargs='+', help='Directories of simulation systems, with _fwd and _bwd fepouts in namd subdirectories')
	ap.add_argument('--prefix', '-', default='dis5A', help='fepout file prefix')
	ap.add_argument('--output', '-o', default='fep_plots.pdf', help='Output file')
	args = ap.parse_args()

	good_dirnames = []
	for dirname in args.dirname:
		if os.path.exists(dirname) is False or os.path.isdir(dirname) is False:
			print('Could not find directory {}.'.format(dirname), file=sys.stderr)
			continue

		fwd_fname = '{}/namd/{}000_fwd.fepout'.format(dirname, args.prefix)
		bwd_fname = '{}/namd/{}000_bwd.fepout'.format(dirname, args.prefix)
		if os.path.exists(fwd_fname) is False:
			print('Could not find {}'.format(fwd_fname))
			continue
		if os.path.exists(bwd_fname) is False:
			print('Could not find {}'.format(bwd_fname))
			continue

		# There should at least be a dis5A000_fwd.fepout
		good_dirnames.append(dirname)

	num_feps = len(good_dirnames)
	print('Processing {} FEPs.'.format(num_feps), file=sys.stderr)

	rows_per_page = 3
	pdf = PdfPages(args.output)

	# Ensure there are enough subplots for all the FEPs
	num_cols = 2
	num_rows = num_feps

	num_pages = int(num_rows / rows_per_page)
	if num_rows % rows_per_page > 0:
		num_pages += 1

	fig, ax = None, None
	# mpl.rcParams['xtick.labelsize'] = 6
	plt.rcParams.update({'font.size': 6})

	for i in range(num_feps):
		page = int(i / rows_per_page)
		offset = i % rows_per_page

		if offset == 0:
			fig, ax = plt.subplots(rows_per_page, num_cols)
			fig.tight_layout(h_pad=2)

		dirname = good_dirnames[i]
		fwd, fwd_dg, lambdas = parse_fepout('{}/namd/{}000_fwd.fepout'.format(dirname, args.prefix))
		bwd, bwd_dg, lambdas = parse_fepout('{}/namd/{}000_bwd.fepout'.format(dirname, args.prefix))
		ax[offset][0].plot(fwd, linewidth=0.2, label='Forward')
		ax[offset][0].plot(bwd, linewidth=0.2, label='Backward')
		ax[offset][0].set_title(dirname, fontsize=8, pad=3)
		# ax[offset][0].set_ylabel('delta-G')
		ax[offset][1].plot(lambdas, fwd_dg, label='Forward')
		ax[offset][1].plot(lambdas, bwd_dg, label='Backward')

		print('{} is on page {}'.format(dirname, page+1))

		if offset == (rows_per_page - 1):
			# Put a legend in the first row
			ax[0][1].legend(loc='best')
			pdf.savefig(fig)

	pdf.close()

if __name__ == '__main__':
	main()
