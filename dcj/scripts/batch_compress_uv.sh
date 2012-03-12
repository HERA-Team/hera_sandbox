#$ -S /bin/bash
ARGS=`pull_args.py $*`
echo compress_uv.py $ARGS
compress_uv.py -x $ARGS

