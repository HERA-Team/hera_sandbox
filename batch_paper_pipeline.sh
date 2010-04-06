#$ -S /bin/bash
ARGS=`pull_args.py $*`
echo paper_pipeline.sh $ARGS
paper_pipeline.sh $ARGS
