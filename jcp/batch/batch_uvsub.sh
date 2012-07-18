#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output/
#$ -N uvadd

ARGS=`pull_args.py $*`
for FILE in ${ARGS}
do 
    uv_addsub_noflags.py ${FILE}F ${FILE}BF
    uv_addsub_noflags.py ${FILE}Fa ${FILE}BBF
    uv_addsub_noflags.py ${FILE}Faa ${FILE}BBB
done
