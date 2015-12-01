#! /bin/bash
# combines multiple outputs from previous lstbinnings
# usage combine_lstbins.sh epoch2/lstbin_July2_v1 epoch3/lstbin_July4_v2
# assumes input directory structure: <input dir>/even/sepX,Y
SEPS=sep0,2
DAYS='even odd'
APPELATION=uvAS
LSTS=`python -c "import numpy as n; hr = 1.0; print ' '.join(['%f_%f' % (d,d+hr/4.) for d in n.arange(0,23.999*hr,hr/4.)])"`
MY_LSTS=`~/src/capo/dcj/scripts/pull_args.py $LSTS`
CALFILE=psa6622_v002
PREFIX='lstbin_epochs2_3_Sep4_v1'
echo creating $PREFIX
mkdir $PREFIX
cd $PREFIX
# 

(
for LST in $MY_LSTS;
do
    #echo working on $LST
    for day in $DAYS;
    do
        echo working on $day
        mkdir $day
        cd $day
        for sep in $SEPS;
        do
            echo working on $sep
            mkdir $sep
            cd $sep
            FILES=
            for LSTDIR in $*;
            do  
                dayfiles=(${LSTDIR}/${day}/${sep}/*${APPELATION})
                shopt -s nullglob
                FILES=$FILES' '${dayfiles[@]}
                shopt -u nullglob
            done

            echo python ~/src/capo/scripts/lstbin_v02.py -a cross -C ${CALFILE}  --lst_res=42.95 --lst_rng=$LST \
    --tfile=600 --altmax=0 --stats=all --median --nsig=6 ${FILES[@]}
             python ~/src/capo/scripts/lstbin_v02.py -a cross -C ${CALFILE}  --lst_res=42.95 --lst_rng=$LST \
    --tfile=600 --altmax=0 --stats=all --median --nsig=6 ${FILES}
            cd ..
        done
        cd ..
    done
done

)



