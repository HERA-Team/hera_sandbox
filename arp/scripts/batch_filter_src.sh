#$ -S /bin/bash
CAL=pgb966_v002
SRCS=J0814+481

ARGS=`pull_args.py $*`
echo filter_src.py -e -C $CAL -s $SRCS -d 3 -r 3 --clean=1e-3 $ARGS
filter_src.py -e -C $CAL -s $SRCS -d 3 -r 3 --clean=1e-3 $ARGS
