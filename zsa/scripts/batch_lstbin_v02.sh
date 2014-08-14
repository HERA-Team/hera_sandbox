#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l paper
#$ -l h_vmem=24G
#LSTS=`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(0,23.999*hr,hr/4.)])"`
LSTS=`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(6,9*hr,hr/4.)])"`
#LSTS=`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(5.25,9.5*hr,hr/4.)])"`
#LSTS=9.250000_9.500000
MY_LSTS=`pull_args.py $LSTS`
CALFILE=psa6240_v003
EXT=$*


echo $MY_LSTS
for LST in $MY_LSTS; do
    echo Working on $LST
    echo /usr/global/paper/capo/scripts/lstbin_v02.py -a cross -C ${CALFILE} -s Sun --lst_res=42.95 --lst_rng=$LST --tfile=600 --altmax=0 --stats=all --median --nsig=3 zen*.${EXT} 
    /usr/global/paper/capo/scripts/lstbin_v02.py -a cross -C ${CALFILE} -s Sun --lst_res=42.95 --lst_rng=$LST --tfile=600 --altmax=0 --stats=all --median --nsig=3 zen*.${EXT} 

done;
