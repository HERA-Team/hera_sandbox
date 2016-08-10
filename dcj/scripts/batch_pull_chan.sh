#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=2G
#$ -j y
#$ -N slice
#$ -o grid_output
#$ -q all.q
PATH=/home/jacobsda/src/anaconda/bin/:$PATH
source activate PAPER
FILES=`~/scripts/pull_args.py $*`
echo pulling data from $FILES
~/scripts/pull_chan.py -p xx -c 100 $FILES
