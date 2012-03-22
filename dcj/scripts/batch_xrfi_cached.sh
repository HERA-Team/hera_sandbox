#! /bin/bash
#$ -N psa332_r
ARGS=`pull_args.py $*`
/data1/paper/arp/scripts/xrfi_simple.py -c 0_143,323,378_387,850_852,900_1023 -t 53 $ARGS
