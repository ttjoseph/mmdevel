#!/usr/bin/env python3
import argparse
import sys
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_fepout(fnames):
	fepenergy, deltas, lambdas = list(), list(), list()
	lines = list()

	if isinstance(fnames, str):
		fnames = [fnames,]

	for fname in fnames:
		print(fname, file=sys.stderr)
		if os.path.exists(fname) is False:
			return None, None, None
		with open(fname) as f:
			lines.extend(f.readlines())

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


# TODO: Return the first matching glob of files
# This will allow us to use (for example) either dis5A000_fwd.fepout, or dis5A???.fepout
def try_globs(*globspecs):
	for globspec in globspecs:
		fnames = glob(globspec)
		if len(fnames) > 0:
			break
	return sorted(fnames)

def main():
	ap = argparse.ArgumentParser(description='Make a big plot of FEP curves')
	ap.add_argument('dirname', nargs='+', help='Directories of simulation systems, with _fwd and _bwd fepouts in namd subdirectories')
	ap.add_argument('--fepdirname', '-f', default='namd', help='Subdirectory under <dirname> containing .fepout files')
	ap.add_argument('--prefix', '-p', default='dis5A', help='fepout file prefix')
	ap.add_argument('--output', '-o', help='Output file')
	args = ap.parse_args()

	good_dirnames = []
	fwd_fnames, bwd_fnames = dict(), dict()
	for dirname in args.dirname:
		if os.path.exists(dirname) is False or os.path.isdir(dirname) is False:
			print('Could not find directory {}.'.format(dirname), file=sys.stderr)
			continue

		fwd_fname = try_globs('{}/{}/{}???_fwd.fepout'.format(dirname, args.fepdirname, args.prefix),
			'{}/{}/{}???.fepout'.format(dirname, args.fepdirname, args.prefix))
		bwd_fname = try_globs('{}/{}/{}???_bwd.fepout'.format(dirname, args.fepdirname, args.prefix))
		if len(fwd_fname) == 0:
			print('Could not find forward fepouts in {}'.format(dirname))
			continue
		else:
			fwd_fnames[dirname] = fwd_fname
		if len(bwd_fname) == 0:
			print('Could not find backward fepouts in {}'.format(dirname))
		else:
			bwd_fnames[dirname] = bwd_fname

		# Only care about this directory if we found a fepout
		good_dirnames.append(dirname)

	num_feps = len(good_dirnames)
	print('Processing {} FEPs.'.format(num_feps), file=sys.stderr)

	if len(good_dirnames) == 0:
		exit(1)

	rows_per_page = 4
	output_filename = args.output or f"{args.prefix}.fep.pdf"
	print(f"Writing plots to: {output_filename}", file=sys.stderr)
	pdf = PdfPages(output_filename)

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
			print('Starting page {}'.format(page+1))
			fig, ax = plt.subplots(rows_per_page, num_cols, figsize=(10, 7.5))
			fig.tight_layout(h_pad=2.0, w_pad=0.2)

		dirname = good_dirnames[i]
		fwd, fwd_dg, lambdas = parse_fepout(fwd_fnames[dirname])
		ax[offset][0].plot(fwd, linewidth=0.2, label='Forward')
		ax[offset][0].set_title(args.prefix, fontsize=8, pad=3)
		ax[offset][1].set_title('{}/{}'.format(os.getcwd(), dirname), fontsize=8, pad=3)
		ax[offset][1].plot(lambdas, np.cumsum(fwd_dg), marker='.', label='Forward: {:.2f} kcal/mol'.format(np.sum(fwd_dg)))
		# ax[offset][1].set_xlabel('λ')
		# ax[offset][1].set_ylabel('ΔG')

		if dirname in bwd_fnames:
			bwd, bwd_dg, lambdas = parse_fepout(bwd_fnames[dirname])
			ax[offset][0].plot(bwd, linewidth=0.2, label='Backward')
			ax[offset][1].plot(lambdas, np.cumsum(bwd_dg), marker='.', label='Backward: {:.2f} kcal/mol'.format(np.sum(bwd_dg)))

		# Show the legend because the dG sums are in it
		ax[offset][1].legend(loc='best')

		if offset == (rows_per_page - 1) or offset == (num_feps - 1):
			pdf.savefig(fig)

	pdf.close()

if __name__ == '__main__':
	main()
