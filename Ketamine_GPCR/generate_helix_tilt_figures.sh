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
SIMULATIONS="5c1m:MuOR/MD_SKP_r2_HSP297 5c1m:MuOR/MD_SKP_r2_HSP297_run2 5c1m:MuOR/MD_SKP_r2_LEU147 $SIMULATIONS"
SIMULATIONS="4djh:KappaOR/MD_SKE_r2 4djh:KappaOR/MD_SKP_r2 $SIMULATIONS"
SIMULATIONS="4djh:KappaOR/MD_SKE_r2_HSP291 4djh:KappaOR/MD_SKP_r2_HSP291 $SIMULATIONS"
SIMULATIONS="4djh:KappaOR/MD_Unbound $SIMULATIONS"
SIMULATIONS="5tvn:5HT2BR/MD_SKE_r2 5tvn:5HT2BR/MD_SKP_r2 $SIMULATIONS"
# SIMULATIONS=
SIMULATIONS="adra2a_1:ADR2A/MD_Unbound adra2a_1:ADR2A/MD_DXM_r1 adra2a_1:ADR2A/MD_DXM_r1_PAL $SIMULATIONS"
SIMULATIONS="adra2a_1:ADR2A/MD_DXM_r1_up2 adra2a_1:ADR2A/MD_AXM adra2a_1:ADR2A/MD_AXM_PAL $SIMULATIONS"

# TM6 helix bend angle specs for angle_over_traj.py, in a fancy bash associative array
declare -A tm6_bend_spec
tm6_bend_spec[5c1m]="PROA:271-292-299"
tm6_bend_spec[4djh]="PROA:268-286-293"
tm6_bend_spec[5tvn]="PROA:321-336-347"
tm6_bend_spec[adra2a_1]="PROA:370-386-397"

# TM3-TM6 distance
declare -A tm3_tm6_dist_spec
tm3_tm6_dist_spec[5c1m]="PROA:169-271"
tm3_tm6_dist_spec[4djh]="PROA:159-266"
tm3_tm6_dist_spec[5tvn]="PROA:157-320"
tm3_tm6_dist_spec[adra2a_1]="PROA:135-370"

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
    #python $MMDEVEL/helix_tilt.py --absolute $helixes ../step5_assembly.xplor_ext.psf ../step5_assembly.namd.pdb prod*.dcd > $gnuplot_script
    # gnuplot $gnuplot_script
    #python $MMDEVEL/csv_stats.py helix-tilt-absolute.csv > $HOME/Dropbox/Ketamine_GPCR_Paper/HelixTilt/${outname}-absolute.csv

    # Forget about non-absolute helix tilts - this is left in for posterity
    #python $MMDEVEL/helix_tilt.py $helixes ../step5_assembly.xplor_ext.psf ../step5_assembly.namd.pdb prod*.dcd > $gnuplot_script
    # gnuplot $gnuplot_script
    #python $MMDEVEL/csv_stats.py helix-tilt.csv > $HOME/Dropbox/Ketamine_GPCR_Paper/HelixTilt/${outname}.csv

    rm -f $gnuplot_script

    # Let's also calculate TM6 helix bend
    pdb=$(echo $sim | cut -d: -f1)
    spec=${tm6_bend_spec[$pdb]}
    python $MMDEVEL/angle_over_traj.py $spec ../step5_assembly.xplor_ext.psf ../step5_assembly.namd.pdb prod*.dcd > tm6-bend.csv
    python $MMDEVEL/csv_stats.py tm6-bend.csv > $HOME/Dropbox/Ketamine_GPCR_Paper/HelixTilt/${outname}-tm6-bend.csv

    # Or maybe the TM3-TM6 distance. Why not? The world is ours
    spec=${tm3_tm6_dist_spec[$pdb]}
    python $MMDEVEL/dist_over_traj.py $spec ../step5_assembly.xplor_ext.psf ../step5_assembly.namd.pdb prod*.dcd > tm3-tm6-dist.csv
    python $MMDEVEL/csv_stats.py tm3-tm6-dist.csv > $HOME/Dropbox/Ketamine_GPCR_Paper/HelixTilt/${outname}-tm3-tm6-dist.csv


done

cd $HERE
