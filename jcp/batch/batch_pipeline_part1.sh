#$ -S /bin/bash
#$ -j y
#$ -o grid_output
#$ -e grid_output
ARGS=`pull_args.py $*`
CAL=psa746_v010
WINDOW='blackman-harris'

for FILE in $ARGS; do   
    #pspec_prep.py -C $CAL --horizon=1.2 --window=none --model --nolstbin ${FILE}
    #pspec_prep.py -C $CAL --horizon=2.0 --window=$WINDOW --model --nolstbin --nophs --nogain ${FILE}B
    pspec_prep.py -C $CAL --nohorizon --window=$WINDOW --model --nolstbin --nophs --nogain --clean=1e-9 ${FILE}BB
done
