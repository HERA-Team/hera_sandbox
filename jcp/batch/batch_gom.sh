#$ -S /bin/bash
#$ -j y
#$-o grid_output
ARGS=`pull_args.py $*`
GOM=1

SRCS1="cas,cyg,crab,vir,Sun"

for FILE in $ARGS; do
    tempgain.py $FILE --gom=$GOM

done
