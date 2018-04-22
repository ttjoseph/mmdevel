#!/bin/bash
#
# Generates the helix tilt plots for the ketamine-GPCR binding affinity project.
# To make the plots look prettier, you will have to modify helix_tilt.py because
# that generates the Gnuplot script.

MMDEVEL=$HOME/mmdevel
B=$HOME

# Format of each entry is <helix-def>:<dir>
# Where dir is immediately under $B and helix-def is a helix-tilt.py helix spec or shortcut
SIMULATIONS="5c1m:MuOR/WithDISU 5c1m:MuOR/MD_SKE_r2.WithDISU 5c1m:MuOR/MD_SKP_r2 5c1m:MuOR/MD_SKP_r2_run2 $SIMULATIONS"
SIMULATIONS="5c1m:MuOR/MD_SKE_r2_HSP297 5c1m:MuOR/MD_SKE_r2_LEU147 $SIMULATIONS"
SIMULATIONS="5c1m:MuOR/MD_SKP_r2_HSP297 5c1m:MuOR/MD_SKP_r2_LEU147 $SIMULATIONS"
SIMULATIONS="4djh:KappaOR/MD_SKE_r2 4djh:KappaOR/MD_SKP_r2 $SIMULATIONS"
SIMULATIONS="4djh:KappaOR/MD_SKE_r2_HSP291 4djh:KappaOR/MD_SKP_r2_HSP291 $SIMULATIONS"
SIMULATIONS="5tvn:5HT2BR/MD_SKE_r2 5tvn:5HT2BR/MD_SKP_r2 $SIMULATIONS"

HERE=$PWD

for sim in $SIMULATIONS
do
    helixes=$(echo $sim | cut -d: -f1)
    dir=$(echo $sim | cut -d: -f2)
    outname=$(echo $dir | sed s#/#-#)
    cd $B/$dir/namd || exit
    echo Working on: $sim prod*.dcd

    # Perhaps the best strategy here is to calculate absolute helix tilts, then calculate differences
    # outside helix_tilt.py. Perhaps in that Excel spreadsheet. Otherwise, we'd have to do unholy
    # merging of structure/coordinates in order to get things to work.
    # We can calculate the helix tilts for 4DKL (and the analogous KOR structure) once, and postprocess.

    gnuplot_script=`mktemp`
    python $MMDEVEL/helix_tilt.py --absolute $helixes ../step5_assembly.xplor_ext.psf ../step5_assembly.namd.pdb prod*.dcd > $gnuplot_script
    # gnuplot $gnuplot_script
    python $MMDEVEL/csv_stats.py helix-tilt-absolute.csv > $HOME/Dropbox/Ketamine_GPCR_Paper/HelixTilt/${outname}-absolute.csv

    #python $MMDEVEL/helix_tilt.py $helixes ../step5_assembly.xplor_ext.psf ../step5_assembly.namd.pdb prod*.dcd > $gnuplot_script
    # gnuplot $gnuplot_script
    #python $MMDEVEL/csv_stats.py helix-tilt.csv > $HOME/Dropbox/Ketamine_GPCR_Paper/HelixTilt/${outname}.csv

    rm -f $gnuplot_script

done

cd $HERE
