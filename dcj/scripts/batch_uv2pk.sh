#$ -S /bin/bash
#$ -V
#$ -j y
#$ -cwd
#$ -o grid_output/
ARGS=`pull_args.py $*`
CAL=psa331_v008_gc
for FILE in $ARGS; do
#    CAL=`python -c "print ['pgb015_v006','pgb015_v007'][int(int('${FILE}'.split('.')[1]) > 2455022)]"`
    echo uv2pk006.py -a cross,-24 -p xx -C $CAL ${FILE}
    uv2pk006.py -a cross,-24 -p xx -C $CAL ${FILE}
done
