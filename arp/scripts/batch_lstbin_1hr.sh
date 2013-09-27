#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l paper
#$ -l h_vmem=4G
LSTS=`python -c "import numpy as n; hr = n.pi/12; print ' '.join(['%f_%f' % (d,d+hr) for d in n.arange(0,24*hr,hr)])"`
MY_LSTS=`pull_args.py $LSTS`

echo $MY_LSTS
for LST in $MY_LSTS; do
    echo Working on $LST
    #lstbin.py -a cross -C psa898_v003 -s Sun --lst_res=42.95 --lst_rng=$LST --tfile=600 $*
    lstbin.py -a cross -C psa898_v003 -s Sun --lst_res=42.95 --lst_rng=$LST --tfile=3600 $*
done;
