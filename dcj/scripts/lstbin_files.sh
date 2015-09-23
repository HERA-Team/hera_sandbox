#! /bin/bash
CAL=psa6622_v001
SUN=n
LSTS=`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(0,23.999*hr,hr/4.)])"`
#LSTBIN=60 #bin size in minutes
#LSTs=`python -c "import numpy; print ' '.join(numpy.arange(0,24,${LSTBIN}/60.).astype('|S8'))"`
for LSTBIN in $LSTS
do
#LST_end=`python -c "print ${LST}+${LSTBIN}*1/60."`
FILES=`~/scripts/lst_select.py --ra=${LSTBIN} --suntime=${SUN} -C ${CAL} $* | wc -l`
echo ${LSTBIN} ${FILES}
done
