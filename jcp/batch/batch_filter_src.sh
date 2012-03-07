#$ -S /bin/bash
#$ -V
#$ -j y
#$ -cwd
#$ -o grid_output/
#$ -N filter_sun
CAL=pgb015_v008
SRCS=Sun

ARGS=`pull_args.py $*`
for FILE in $ARGS; do
    filter_src.py -C $CAL -s $SRCS -d 3 --clean=1e-4 $FILE
done
