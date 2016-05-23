#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -e grid_output
#$ -o grid_output
#$ -l paper
#$ -l h_vmem=16G

FILES=`pull_args.py $*`
CALFILE=hsa7458_v000

for f in ${FILES}; do
    echo omni_run_hera.py -p xx -C hsa7458_v000 --calpar=/data4/paper/HERA2015/fcgains.xx.npz --omnipath=/data4/paper/HERA2015/2457458/omnifiles/ --ba=81 ${f}
    omni_run_hera.py -p xx -C hsa7458_v000 --calpar=/data4/paper/HERA2015/fcgains.xx.npz --omnipath=/data4/paper/HERA2015/2457458/omnifiles/ --ba=81 ${f}
done;
