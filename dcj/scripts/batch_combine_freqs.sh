#$ -S /bin/bash
#$ -V
#$ -j y
#$ -cwd
#$ -o grid_output/
#$ -N avg_freqs

FILES=`pull_args.py $*`
combine_freqs.py -n 512 ${FILES}
