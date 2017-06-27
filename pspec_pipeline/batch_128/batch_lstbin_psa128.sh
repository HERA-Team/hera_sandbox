#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -j y
#$ -l h_vmem=18G
#$ -N LSTBIN
export PATH=/home/saulkohn/src/anaconda/bin:$PATH
source activate PAPER
LSTS=`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(0,23.999*hr,hr/4.)])"`
#LSTS=`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(6,9*hr,hr/4.)])"`
#LSTS=`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(5.25,9.5*hr,hr/4.)])"`
#LSTS=9.250000_9.500000
MY_LSTS=`pull_args.py $LSTS`
CALFILE=psa6622_v002
PREFIX=lstbin_Apr12_V_fgsub
PATH2CAPO=/home/saulkohn/ReposForCanopy/capo
echo $MY_LSTS
echo mkdir ${PREFIX}
mkdir ${PREFIX}
cd ${PREFIX}
for LST in $MY_LSTS; do
    echo Working on $LST
    echo working on even files
    mkdir even
    cd even
    echo ${PATH2CAPO}/scripts/lstbin_v02.py -a cross -C ${CALFILE} --lst_res=42.95 --lst_rng=$LST \
    --tfile=600 --altmax=0 --stats=all --median --nsig=3 `python ${PATH2CAPO}/scripts/select_file_parity.py $*`
    python ${PATH2CAPO}/scripts/lstbin_v02.py -a cross -C ${CALFILE}  --lst_res=42.95 --lst_rng=$LST \
    --tfile=600 --altmax=0 --stats=all --median --nsig=3 `python ${PATH2CAPO}/scripts/select_file_parity.py $*`
    cd ..

    echo working on odd files
    mkdir odd
    cd odd
    echo ${PATH2CAPO}/scripts/lstbin_v02.py -a cross -C ${CALFILE} --lst_res=42.95 --lst_rng=$LST \
    --tfile=600 --altmax=0 --stats=all --median --nsig=3 `python ~/scripts/select_file_parity.py $*`
    python ${PATH2CAPO}/scripts/lstbin_v02.py -a cross -C ${CALFILE} --lst_res=42.95 --lst_rng=$LST \
    --tfile=600 --altmax=0 --stats=all --median --nsig=3 `python ${PATH2CAPO}/scripts/select_file_parity.py --odd $*`
    cd ..

done;
