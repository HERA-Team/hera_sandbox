#$ -S /bin/bash
#CAL=pgb015_v005
CAT=test_cat_gal_cut_5.a
CAL=pgb015_v005_gc
#SRCS=`srclist.py -s 30/.150 -x 800/.150 --dec=-10_90`
#SRCS=`python -c "import $CAL ; print ' '.join(${CAL}.src_prms.keys())"`
#SRCLIST="216 286 239"
SRCLIST=`srclist.py -C ${CAL} -s 10/0.15 --cat=${CAT} --dec=10_60`
#SRCS="cyg cas crab vir Sun"
#SRCLIST=`python -c "print '$SRCS'.replace(' ',',')"`
#CAT=paper_deep
#CAT=three_cr
#SRCLIST=`srclist.py -s all --cat=${CAT} `
SRCS=`pull_args.py ${SRCLIST}`
#FILES=`lst_select.py -s ${SRCS} --cat=${CAT} $*`
DLYWID=10
DRTWID=10
MINUV=1
CLN=5e-4
#MINUV=10
echo SRCS=${SRCS}
echo $*
for SRC in $SRCS ; do
    mkdir $SRC
    cd $SRC
    cp $* 
    #echo mdlvis.py -C $CAL -s $SRCLIST -m sub
    #mdlvis.py -C $CAL -s $SRCLIST -m sub $ARG
    mdlvis.py -s `srclist.py -s ${SRCS} --exclude=${SRC}` $* -C ${CAL} -f 
    SRCFILES=`lst_select.py -C ${CAL} -s ${SRC} --cat=${CAT} $*`
    echo Processing $SRC,${SRCFILES}
    filter_src.py -e -C $CAL -s $SRC --cat=${CAT} -d $DLYWID -r $DRTWID --clean=$CLN ${SRCFILES}
    FSRCS=`python -c "print ' '.join([s+'.e_${SRC}' for s in '${SRCFILES}'.split()])"`
        #mdlvis.py -C $CAL -s $SRC -m add ${ARG}s.e_${SRC}
    echo $FSRCS
    beamform.py -C $CAL -p yy --minuv=$MINUV -s $SRC $FSRCS --cat=${CAT}
    FBSRCS=`python -c "print ' '.join([s+'*.bm_${SRC}' for s in '${FSRCS}'.split()])"`
        #combine_freqs.py -n 16 -u $FBSRCS
    sum_all_integrations.py  ${FBSRCS}
    FBASRCS=`python -c "print ' '.join([s + '*A' for s in '${FBSRCS}'.split()])"`
    combine_freqs.py -d -u -n 32 ${FBASRCS}
    ~/scripts/plot_uv.py  -pyy -t0 -mlin --chan_axis=physical *bm_${SRC}Am | tee ${CAL}_${SRC}.txt
done
