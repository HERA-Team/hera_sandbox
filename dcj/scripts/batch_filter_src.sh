#$ -S /bin/bash
#$ -V
#$ -j y
#$ -cwd
#$ -o grid_output/
#$ -N filter_sun
CAL=psa331_v008_gc
SRCS=Sun

ARGS=`pull_args.py $*`
echo filter_src.py -e -C $CAL -s $SRCS -d 5 -r 5 --clean=5e-4 $ARGS
filter_src.py -e -C $CAL -s $SRCS -d 5 -r 5 --clean=5e-4 $ARGS
