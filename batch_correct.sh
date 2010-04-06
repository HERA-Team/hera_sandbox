#$ -S /bin/bash
ARGS=`pull_args.py $*`
CORRECT=correct_psa113.py
BPSCALE=275560000

for FILE in $ARGS; do
    $CORRECT $FILE
    apply_bp.py -s $BPSCALE ${FILE}c
    sum_integrations.py -n 10 ${FILE}cb
done
