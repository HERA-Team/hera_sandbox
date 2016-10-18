#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -e grid_output
#$ -o grid_output
#$ -l paper
#$ -l h_vmem=5G


FILES=`pull_args.py $*`
CALFILE=hsa7458_v000

for f in ${FILES}; do
    echo pull_chan.py -p xx -c 40_60,110_150 --bls=auto ${f}
    pull_chan.py -p xx -c 40_60,110_150 --bls=auto ${f}
    mv ${f}.npz autochans/.
done




