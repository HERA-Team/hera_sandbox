#$ -S /bin/bash
#CAL=pgb015_v004
CAL=pgb015_v006
SPECSRCS="cyg cas crab vir"
#ALLSRCS=`srclist.py -s 20/.150 -x 500/.150 --dec=8_68 --ra=10_5`
#ALLSRCS="cyg cas crab vir"
#ALLSRCS="${ALLSRCS} her hyd cyg cas crab vir"
ALLSRCS=`python -c "import $CAL ; print ' '.join(${CAL}.src_prms.keys())"`
#ALLSRCS="J0505+381 J0543+495 J0655+541 J0814+481 J0921+454 J1524+543 J1535+554 J1659+470 J1724+506 J1830+484"
SRCS=`pull_args.py $ALLSRCS`
#CH=240_720
CH=60_180
#ANTS=cross
#ANTS=cross,-0_5,-0_6,-1_8,-3_11,-3_12,-4_15,-5_6,-12_13
ANTS=cross,-0_5,-0_6,-3_12,-5_6

for SRC in $SRCS ; do
    echo $SRC
    #if [[ $SPECSRCS == *$SRC* ]] ; then
        meas_src.py -p yy -C $CAL -c $CH -a $ANTS --altmin=30 -s $SRC -o ${CAL}_srcspec_${SRC}.png *.e_${SRC}s
    #else
        #meas_src.py -p yy -C $CAL -c $CH -a $ANTS --altmin=30 -s $SRC -o ${CAL}_srcspec_${SRC}.png *.uvdddddrm
    #fi
done
