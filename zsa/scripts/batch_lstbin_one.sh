#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l paper
#$ -l h_vmem=12G
#LSTS=`python -c "import numpy as n; hr = n.pi/12; print ' '.join(['%f_%f' % (d,d+hr) for d in n.arange(7*hr,9*hr,hr)])"`
#LST=1.832596_2.094395 
#LST=2.094395_2.356194 
LST=2.356194_2.617994
ITERS=`python -c "for i in range(6): print i"`
ITERS=`pull_args.py $ITERS`
FILES=(*.uvcRREcAzxCPBxR)
CALFILE=psa6240_v003

#564 files broken up into 94 each. 6 blocks.
for i in $ITERS; do
    st=0
    cnt=$[ $[${i}+1]*94 ]
    ARR=( "${FILES[@]:${st}:${cnt}}" )
    echo lstbin.py -a cross -C ${CALFILE} -s Sun --lst_res=42.95 --lst_rng=$LST --tfile=3600 `echo ${ARR[@]}`
    lstbin.py -a cross -C ${CALFILE} -s Sun --lst_res=42.95 --lst_rng=$LST --tfile=3600 --ending=${i} `echo ${ARR[@]}`
done;
