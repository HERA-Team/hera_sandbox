#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l paper
#$ -j y
#$ -l h_vmem=12G
#$ -N LSTBIN

#export PATH=/home/jacobsda/src/anaconda/bin:$PATH
#source activate PAPER

LSTS=`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(0,23.999*hr,hr/4.)])"` #96 lst bins from 0-24hr
#LSTS=`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(6,9*hr,hr/4.)])"`
#LSTS=`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(5.25,9.5*hr,hr/4.)])"`
#LSTS=9.250000_9.500000

MY_LSTS=`pull_args.py $LSTS`
CALFILE=psa6622_v003 #Aaron's fast calfile
PREFIX=lstbin_full


echo $MY_LSTS
echo mkdir ${PREFIX}
mkdir ${PREFIX}
cd ${PREFIX}

for LST in $MY_LSTS; do
    echo Working on $LST
    echo working on even files
    mkdir even
    cd even
    #echo ~/capo/scripts/lstbin_v02.py -a cross -C ${CALFILE} --lst_res=31.65 --lst_rng=$LST \
    #--tfile=600 --altmax=0 --stats=all --median --nsig=3 -s Sun `python ~/capo/ctc/scripts/select_file_parity.py $*`
    python ~/capo/scripts/lstbin_v02.py -a cross -C ${CALFILE}  --lst_res=31.65 --lst_rng=$LST \
    --tfile=600 --altmax=0 --stats=all --median --nsig=3 -s Sun `python ~/capo/ctc/scripts/select_file_parity.py $*`
    #python ~/capo/scripts/lstbin_v02.py -a cross -C ${CALFILE}  --lst_res=31.65 --lst_rng=$LST \
    #--tfile=600 --altmax=0 --stats=all --median --nsig=3 -s Sun `python ~/capo/sak/scripts/alternate_djd_select.py $*`
    cd ..
    echo Working on odd files
    mkdir odd
    cd odd
    #echo ~/capo/scripts/lstbin_v02.py -a cross -C ${CALFILE} -s Sun --lst_res=31.65 --lst_rng=$LST \
    #--tfile=600 --altmax=0 --stats=all --median --nsig=3 -s Sun `python ~/capo/ctc/scripts/select_file_parity.py $*`
    python ~/capo/scripts/lstbin_v02.py -a cross -C ${CALFILE} --lst_res=31.65 --lst_rng=$LST \
    --tfile=600 --altmax=0 --stats=all --median --nsig=3 -s Sun `python ~/capo/ctc/scripts/select_file_parity.py --odd $*`
    #python ~/capo/scripts/lstbin_v02.py -a cross -C ${CALFILE}  --lst_res=31.65 --lst_rng=$LST \
    #--tfile=600 --altmax=0 --stats=all --median --nsig=3 -s Sun `python ~/capo/sak/scripts/alternate_djd_select.py --odd $*`
    cd ..

done;
