#$ -S /bin/bash
CAL=psa112_v005
ARGS=`pull_args.py $*`
SRCS=Sun,pic,hyd,for,J2214-170,cyg,crab,J1615-605,1932-464,J2154-692,J2358-605,cen
CATS=helm,misc,culgoora
#SRCS=${SRCS},`srclist.py -s 10/.150 -x 500/.150 --dec=20_50 --divstr=,`
#SRCS=`python -c "import $CAL ; print ','.join(${CAL}.src_prms.keys())"`

echo mdlvis.py -C $CAL -s $SRCS -f $ARGS
mdlvis.py -C $CAL -s $SRCS --cat=$CATS -f $ARGS

