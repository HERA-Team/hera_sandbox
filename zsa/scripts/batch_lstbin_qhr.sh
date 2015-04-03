#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l paper
#$ -l h_vmem=12G
LSTS=`python -c "import numpy as n; hr = n.pi/12; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(0,24*hr,hr/4.)])"`
MY_LSTS=`pull_args.py $LSTS`
CALFILE=psa6240_v003

echo $MY_LSTS
for LST in $MY_LSTS; do
    echo Working on $LST
    #lstbin.py -a cross -C psa898_v003 -s Sun --lst_res=42.95 --lst_rng=$LST --tfile=600 $*
    #lstbin.py -a cross -C psa898_v003 -s Sun --lst_res=42.95 --lst_rng=$LST --tfile=3600 $*
    #lstbin.py -a cross -C psa898_v003 -s Sun --lst_res=42.95 --lst_rng=$LST --tfile=3600 *EzCPSBxR
    lstbin.py -a cross -C ${CALFILE} -s Sun --lst_res=42.95 --lst_rng=$LST --tfile=3600 *.uvcRREcACP
done;
