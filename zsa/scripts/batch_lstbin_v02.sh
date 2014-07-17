#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l paper
#$ -l h_vmem=12G
LSTS=`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr) for d in n.arange(0,24*hr,hr)])"`
MY_LSTS=`pull_args.py $LSTS`
CALFILE=psa6240_v003


echo $MY_LSTS
for LST in $MY_LSTS; do
    echo Working on $LST
    echo lstbin_v02.py -a cross -C ${CALFILE} -s Sun --lst_res=42.95 --lst_rng=$LST --tfile=600 --altmax=0 --stats=all --median --nsig=3 *.uvcRREcAzxCPB
    lstbin_v02.py -a cross -C ${CALFILE} -s Sun --lst_res=42.95 --lst_rng=$LST --tfile=600 --altmax=0 --stats=all --median --nsig=3 *.uvcRREcAzxCPB
done;
