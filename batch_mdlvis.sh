#$ -S /bin/bash
CAL=psa112_v003
ARGS=`pull_args.py $*`
#SRCS=${SRCS},`srclist.py -s 10/.150 -x 500/.150 --dec=20_50 --divstr=,`
#SRCS=`python -c "import $CAL ; print ','.join(${CAL}.src_prms.keys())"`
SRCS=cyg,crab,hyd,pic,for,J1615-605,J1935-461,J2154-692

echo mdlvis.py -C $CAL -s $SRCS $ARGS
mdlvis.py -C $CAL -s $SRCS $ARGS

