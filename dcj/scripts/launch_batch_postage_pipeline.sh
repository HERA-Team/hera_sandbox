#! /bin/bash
. ~/.bashrc
cd /data1/paper/gb/
UVCBTS=zen.245487[1,2].*.uvcbt
LSTS_ALL=`lst_select.py -C pgb966_v003 --ra=10_14 -d - $UVCBTS`
#NDAYS=`echo $LSTS_ALL |python -c "import sys; print len(sys.stdin.readlines()[0].strip().split('-'))-1"`
NDAYS=`echo "$LSTS_ALL" | sed -n 's/\(-\)/\1/p' | awk 'END {print NR}'`
echo "starting qsub with LSTS_ALL="$LSTS_ALL
qsub  -V -cwd -t 1:$NDAYS:1 /data1/paper/jacobsda/scripts/batch_postage_pipeline.sh $UVCBTS