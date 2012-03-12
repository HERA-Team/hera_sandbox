#! /bin/bash
. ~/.bashrc
cd /data1/paper/gb/
#UVCBS=zen.245487[1,2].*.uvcb
UVCBS=`/data1/paper/jacobsda/scripts/find_difference.py '*.e_Sunm' '*.uvcb'`
LSTS_ALL=`lst_select.py -C pgb966_v003 --sun -d - $UVCBS`
#NDAYS=`echo $LSTS_ALL |python -c "import sys; print len(sys.stdin.readlines()[0].strip().split('-'))-1"`
NDAYS=`echo "$LSTS_ALL" | sed -n 's/\(-\)/\1/p' | awk 'END {print NR}'`
echo "starting qsub with LSTS_ALL="$LSTS_ALL
qsub  -V -cwd -t 1:$NDAYS:1 /data1/paper/jacobsda/scripts/batch_sun_postage_pipeline.sh $UVCBS
