#$ -S /bin/bash
#$ -V
#$ -j y
#$ -cwd
#$ -o grid_output/
#$ -N uvadd

ARGS=`pull_args.py $*`
for FILE in ${ARGS}
do 
    uv_addsub_noflags.py ${FILE}F ${FILE}BF
    uv_addsub_noflags.py ${FILE}Fa ${FILE}BB
done