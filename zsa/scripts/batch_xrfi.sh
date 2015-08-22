#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=4G
#$ -l paper
ARGS=`pull_args.py $*`

for FILE in $ARGS; do     
    echo xrfi_simple.py -n 3 ${FILE} 
    xrfi_simple.py -n 3 ${FILE} 
done
