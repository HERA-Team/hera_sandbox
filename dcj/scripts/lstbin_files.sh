#! /bin/bash

CAL=psa746_v011
SUN=n
LSTBIN=60 #bin size in minutes
LSTs=`python -c "import numpy; print ' '.join(numpy.arange(0,24,${LSTBIN}/60.).astype('|S8'))"`
for LST in $LSTs
do
LST_end=`python -c "print ${LST}+${LSTBIN}*1/60."`
FILES=`~/scripts/lst_select.py --ra=${LST}_${LST_end} --suntime=${SUN} -C ${CAL} $* | wc -l`
echo ${LST} ${FILES}
done
