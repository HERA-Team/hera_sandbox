#$ -S /bin/bash
ARGS=`pull_args.py $*`
echo apply_bp.py $ARGS
apply_bp.py $ARGS
