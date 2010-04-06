#$ -S /bin/bash
ARGS=`pull_args.py $*`
echo precal_pipeline.sh $ARGS
precal_pipeline.sh `pull_args.py $*`
