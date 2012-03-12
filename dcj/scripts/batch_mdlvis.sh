#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output/
#$ -j y
CAL=psa331_v005_gc
ARGS=`pull_args.py $*`
#SRCS=Sun,pic,hyd,for,J2214-170,cyg,crab,J1615-605,1932-464,J2154-692,J2358-605,cen
#CATS=helm,misc,culgoorai
#CATS=pgb015_paper_r_155_good_arp_excluded
CATS=sumss_top1000
SRCS=`srclist.py -s 10/0.15 --ra=9_21 --cat=${CATS} --divstr=, -C ${CAL}`
#SRCS=`python -c "import $CAL ; print ','.join(${CAL}.src_prms.keys())"`

echo mdlvis.py -C $CAL -s $SRCS -f --cat=${CATS} $ARGS
/usr/global/paper/bin/mdlvis.py -C $CAL -s $SRCS --cat=$CATS -f $ARGS
#./mdlvis_time_travel.py -C $CAL -s $SRCS --cat=$CATS -f $ARGS --dt=2
