#$ -S /bin/bash
ARGS=`pull_args.py $*`
echo xrfi_simple.py -c 725_1023 -n 3.5 --combine -t 10 $ARGS
xrfi_simple.py -c 725_1023 -n 3.5 --combine -t 10 $ARGS
