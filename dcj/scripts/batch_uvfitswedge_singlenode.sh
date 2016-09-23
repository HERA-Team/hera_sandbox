#! /bin/bash
#run uvfitswedge.py on a bunch of data
NTASKS=10
PIDS=
Is=`python -c "print ' '.join(map(str,range(1,$NTASKS+1)))"`
for I in $Is;
do
MYDATA=`~/scripts/pull_args.py -t 1:${NTASKS} --taskid=${I} $*`
echo 'node'${I}' working on: '${MYDATA}
python ~/scripts/wedge.py ${MYDATA} | tee wedge_${I}.log &
PIDS="${PIDS} $!"
done
echo "done launching ${NTASKS} jobs, waiting on pids"
echo $PIDS
wait $PIDS
