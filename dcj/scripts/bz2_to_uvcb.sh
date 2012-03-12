#$ -S /bin/bash


month=$(date +%m)
year=$(date +%y)

cd /data1/paper/gb/${year}-${month}
TARS=`/data1/paper/jacobsda/scripts/find_difference.py '*.uvcb' '*.bz2'`
echo 'Decompressing: '
echo $TARS
compress_uv.py -x $TARS

UVS=`/data1/paper/jacobsda/scripts/find_difference.py '*.uvcb' '*.uv'`
echo 'Correcting:'
echo $UVS
python ../arp/scripts/correct_pgb966.py $UVS
rm -r $UVS

UVCS=`/data1/paper/jacobsda/scripts/find_difference.py '*.uvcb' '*.uvc'`
echo 'Applying bp:'
echo $UVCS
apply_bp.py $UVCS
rm -r $UVCS

#UVCBS=`/data1/paper/jacobsda/scripts/find_difference.py '*.uvcbt' '*.uvcb'`
#echo 'Applying bp:'
#echo $UVCBS
#tempgain.py $UVCS
#rm -r $UVCBS

