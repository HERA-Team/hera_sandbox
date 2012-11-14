#$ -S /bin/bash
#$ -j y
#$ -o grid_output
#$ -e grid_output
ARGS=`pull_args.py $*`
CAL=psa746_v010
WINDOW='blackman-harris'

for FILE in $ARGS; do   
    pspec_to_npz_nored.py -C $CAL -p xx ${FILE}
done
