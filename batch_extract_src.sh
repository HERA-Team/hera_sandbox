#$ -S /bin/bash
CAL=pgb015_v005
#SRCS=`srclist.py -s 30/.150 -x 800/.150 --dec=-10_90`
SRCS=`python -c "import $CAL ; print ' '.join(${CAL}.src_prms.keys())"`
#SRCS="cyg cas crab vir"
ARGS=`pull_args.py $*`
DLYWID=5
DRTWID=5
MINUV=20
#MINUV=10

for ARG in $ARGS ; do
    for SRC in $SRCS ; do
        echo Processing $ARG, $SRC
        filter_src.py -e -C $CAL -s $SRC -d $DLYWID -r $DRTWID --clean=1e-3 $ARG
        mdlvis.py -C $CAL -s $SRC -m add ${ARG}.e_${SRC}
        beamform.py -C $CAL -p yy --minuv=$MINUV -s $SRC $ARG.e_${SRC}s
        #combine_freqs.py -d -u -n 32 $ARG.e_${SRC}s.bm_$SRC
    done
done
