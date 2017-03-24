#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=6G
#$ -l scratch=3G
#$ -l paper
#$ -N OMNI_RUN_BANDPASS
ARGS=`pull_args.py $*`

echo ${ARGS}

for f in ${ARGS}; do
    echo working on ${fi}...
    ~/capo/omni/omni_run_bandpass.py ${f} --save
done

