import sys
import argparse
from tempfile import NamedTemporaryFile
import MDAnalysis as mda
from MDAnalysis.analysis import align
from pmda.rms import RMSF

def main():
	ap = argparse.ArgumentParser(description='Calculate RMSF of protein, by residue')
	ap.add_argument('psf', help='Topology file')
	ap.add_argument('dcd', nargs='+', help='Trajectory file(s)')
	ap.add_argument('--skip-frames', type=int, default=0, help='Number of frames to skip at start of trajectory')
	args = ap.parse_args()

	print('Loading system. I hope you have a lot of RAM...', file=sys.stderr)
	u = mda.Universe(args.psf, args.dcd, in_memory=True)

	print("""Aligning the trajectory first, into a temporary file on disk.
This is due to an implementation quirk of the MDAnalysis parallel module.
NOTE: If the protein is split across periodic cells then this will all be wrong!""", file=sys.stderr)

	protein = u.select_atoms("protein")
	# Fit to the initial frame to get a better average structure
	# (the trajectory is changed in memory)
	prealigner = align.AlignTraj(u, u, select="protein and name CA",
	                             in_memory=True).run()

	# ref = average structure
	ref_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
	# Make a reference structure (need to reshape into a
	# 1-frame "trajectory").
	ref = mda.Merge(protein).load_new(ref_coordinates[:, None, :],
	                                  order="afc")

	aligner = align.AlignTraj(u, ref, select="protein and name CA",
                          in_memory=True).run()

	with NamedTemporaryFile(suffix='.dcd') as tmpdcd:
		# need to write the trajectory to disk for PMDA 0.3.0 (see issue #15)
		with mda.Writer(tmpdcd.name, n_atoms=u.atoms.n_atoms) as W:
		    # Skip frames by only writing out non-skipped frames
		    for ts in u.trajectory[args.skip_frames:]:
	        	W.write(u.atoms)		
		u = mda.Universe(args.psf, tmpdcd.name)
		calphas = protein.select_atoms("protein and name CA")

		rmsfer = RMSF(calphas).run(n_blocks=8)

		data = zip(calphas.resnums, rmsfer.rmsf)

		print('resnum,rmsf')
		for d in data:
			print(f'{d[0]},{d[1]}')

if __name__ == '__main__':
	main()