#!/bin/bash
# Conveniently submit a Gaussian job to the friendly local SLURM.
# Usage: submit-g09.sh whatever.gau
nproc=`grep -i nproc $1 | cut -d= -f2`
sbatch <<FOO
#!/bin/bash
#SBATCH -J $1
#SBATCH -n $nproc

g09 $1
echo "Done with $1 in \$PWD" | ~/bin/slack-ttjoseph
FOO
