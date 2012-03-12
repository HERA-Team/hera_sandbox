#$ -S /bin/bash
#$ -N psa332_e
#$ -j y
#$ -o grid_output

#CAL=pgb015_v005

CAL=psa331_v005_gc
#SRCS=`srclist.py -s 30/.150 -x 800/.150 --dec=-10_90`
#SRCS=`python -c "import $CAL ; print ' '.join(${CAL}.src_prms.keys())"`
#SRCS="hyd vir her sgr"
SRCS=Sun
#SRCLIST=`python -c "print '$SRCS'.replace(' ',',')"`
ARGS=`pull_args.py $*`
DLYWID=5
DRTWID=5
MINUV=10
CLN=5e-4
#MINUV=10

for ARG in $ARGS ; do
#    echo mdlvis.py -C $CAL -s $SRCLIST -m sub
#    mdlvis.py -C $CAL -s $SRCLIST -m sub $ARG
    for SRC in $SRCS ; do
        echo Processing $ARG, $SRC
        filter_src.py -e -C $CAL -s $SRC -d $DLYWID -r $DRTWID --clean=$CLN ${ARG}
        uv_addsub.py --sub $ARG ${ARG}.e_${SRC}
#        mdlvis.py -C $CAL -s $SRC -m add ${ARG}s.e_${SRC}
#        uv_addsub.py --sub $ARG
        #beamform.py -C $CAL -p yy --minuv=$MINUV -s $SRC $ARG.e_${SRC}s
        ##combine_freqs.py -d -u -n 32 $ARG.e_${SRC}s.bm_$SRC
    done
done
ssh djm -n Shredder -m "$JOBID is finished"