#!/bin/bash

# Walk an lstbinned dataset and generate a noise realization 
# usage:
# batch_noise_sim.sh  <input lstbinned data>  <name of simulation to create>
# batch_noise_sim.sh lstbin_data lstbin_noise_realization


seps='sep0,1 sep1,1 sep-1,1'
cal=psa6240_v003
indir=$1
outdir=$2
#noiselevel=3.85 #Jy
Trcvr=200
appelation='uvGA'
chan='101'

scriptsdir=/home/djacobs/scripts


declare -A ants
ants[sep0,1]=0_44
ants[sep1,1]=1_3
ants[sep-1,1]=1_48

days='even odd'

printf 'This is Batch Noise Sim using:\n'
echo 'Trcvr = ' ${Trcvr} 
echo 'operating on lstbinned data:' $indir

sleep 1.5

mkdir $outdir
for path in $days; do
    mkdir $outdir/$path
    for sep in $seps; do
        mkdir $outdir/$path/$sep
        files=$(ls -d $indir/$path/$sep/*.${appelation})
        printf 'generating noise'
        echo  mdlvis.py --Trcvr=${Trcvr} --inttime=43 $files --outdir=${outdir}/$path/$sep
        "${scriptsdir}/sim_noise.py" --Trcvr=${Trcvr} --inttime=43 $files --outdir=${outdir}/$path/$sep
    done
done
