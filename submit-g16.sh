#!/bin/bash
# Conveniently submit a Gaussian job to the friendly local SLURM.
# Usage: submit-g16.sh whatever.gau
# nproc=`grep -i nproc $1 | cut -d= -f2`
sbatch <<FOO
#!/bin/bash
#SBATCH -J $1
#SBATCH -n 32

export g16root=/nas1/sw
export PGI_FASTMATH_CPU=haswell
. \$g16root/g16/bsd/g16.profile

g16 $1
echo "Done with $1 in \$PWD" | /nas1/sw/bin/slack-ttjoseph
FOO
