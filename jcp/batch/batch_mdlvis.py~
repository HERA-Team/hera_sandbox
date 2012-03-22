#$ -S /bin/bash
#$ -V
#$ -j y
#$ -cwd
#$ -o grid_output/
#$ -N sub_Sun
CAL=pgb322_v008_gc

ARGS=`pull_args.py $*`
for FILE in ${ARGS}
do 
    echo uv_addsub.py --sub ${FILE} ${FILE}.e_Sun
    uv_addsub.py --sub ${FILE} ${FILE}.e_Sun
done