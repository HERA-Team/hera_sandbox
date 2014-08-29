#! /bin/bash

DATADIR=/home/obs/data
#FOLIO=damo@folio.sas.upenn.edu
FOLIO=djacobs@enterprise.sese.asu.edu
#FolioBase=/data4/raw_data/EoR2013
FolioBase=/data3/PAPER/psa128/
tjd0=''
today=`date -u "+%F"`

python /home/obs/Share/redux/dailyRFIreport.py ${DATADIR}/zen.*yy*.npz -o /home/obs/${today}.png
scp /home/obs/${today}.png damo@enterprise.sese.asu.edu:/data2/PAPER/psa128/RFI/ 

for FILE in `ls -1d $DATADIR/zen*`
do
    #Read filename to get truncated julian date.
    tjd=${FILE##*zen.245}
    tjd=${tjd%%.*}

    #start making directory structure
    oddoreven=$(( $tjd % 2 ))
    #for now, send everything to pot0
    oddoreven=0
    #if directories don't exist, make them. Use ssh as little as possible...
    if [[ "$tjd" -ne "$tjd0" ]]
    then
        tjd0=$tjd
        echo New directory psa${tjd}
        ssh pot${oddoreven} test -e /data${oddoreven}/psa${tjd} \
            || ssh -q -o ConnectTimeout=3 pot${oddoreven} "mkdir /data${oddoreven}/psa${tjd}"
        ssh $FOLIO test -e $FolioBase/psa${tjd} \
            || ssh -q -o ConnectTimeout=3 $FOLIO "mkdir ${FolioBase}/psa${tjd}"
    fi
    echo rsync -aruP $FILE pot${oddoreven}:/data${oddoreven}/psa${tjd}/
    rsync -aruP $FILE pot${oddoreven}:/data${oddoreven}/psa${tjd}/
    if [[ $FILE == *zen.24*E* ]]
    then
        echo rsync -aruP -e "ssh -c arcfour128" $FILE ${FOLIO}:${FolioBase}/psa${tjd}/
        rsync -aruP -e "ssh -c arcfour128" $FILE ${FOLIO}:${FolioBase}/psa${tjd}/
    fi
    rm -r $FILE
done
