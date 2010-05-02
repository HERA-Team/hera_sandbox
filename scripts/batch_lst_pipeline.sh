#$ -S /bin/bash
ARGS=`pull_args.py $*`
echo lst_pipeline.sh $ARGS
lst_pipeline.sh $ARGS
