#$ -S /bin/bash
#$ -j y
#$ -N bz2_to_uvcbt
#$ -o /data1/paper/gb/nightly_pipeline_logs/
cd /data1/paper/gb
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

UVCBS=`/data1/paper/jacobsda/scripts/find_difference.py '*.uvcbt' '*.uvcb'`
echo 'Applying bp:'
echo $UVCBS
tempgain.py $UVCBS
rm -r $UVCBS

