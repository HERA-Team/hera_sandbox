#$ -S /bin/bash
ARGS=`pull_args.py $*`

for FILE in $ARGS; do
    #/data1/paper/arp/scripts/xrfi_simple.py --dt=3.5 --combine -t 80 ${FILE}
    /data3/paper/pober/capo/jcp/scripts/xrfi_simple.py -c 0_130,755_777,1540,1704,1827,1868,1901_2047 --dt=4 --df=6
done
