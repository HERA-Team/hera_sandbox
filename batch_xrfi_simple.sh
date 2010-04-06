#$ -S /bin/bash
ARGS=`pull_args.py $*`
echo xrfi_simple.py --dt=2 --df=3 -c 0_59,805_885,925_1023 --combine -t 20 $ARGS
xrfi_simple.py --dt=2 --df=3 -c 0_59,805_885,925_1023 --combine -t 20 $ARGS
