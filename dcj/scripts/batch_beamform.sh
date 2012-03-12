#$ -S /bin/bash
CAL=psa112_v005
MINUV=20
#SRCS=`python -c "import $CAL ; print ' '.join(${CAL}.src_prms.keys())"`
SRCS = `srclist.py -c all --cat=paper`
SRCS=`pull_args.py $SRCS`

for SRC in $SRCS ; do
    echo beamform.py -C $CAL -p yy --minuv=$MINUV -s $SRC 
    beamform.py -C $CAL -p yy --minuv=$MINUV -s $SRC *.e_${SRC}s
    echo combine_freqs.py -d -u -n 16 *.e_${SRC}s.bm_${SRC}
    combine_freqs.py -d -u -n 16 *.e_${SRC}s.bm_${SRC}
done
