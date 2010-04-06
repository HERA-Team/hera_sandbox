#$ -S /bin/bash
#---------------------------------------------------------------------------
# Admittedly complicated header to first split files between scripts, and if
# more nodes remain, split data within a file across nodes using decimation
# parameters.
ARGS=`pull_args.py -w $*`
NREPEAT=$(((${SGE_TASK_LAST}-${SGE_TASK_FIRST}+1)/$#))
NREMAIN=$(((${SGE_TASK_LAST}-${SGE_TASK_FIRST}+1)%$#))
if [[ $((${SGE_TASK_ID}%$#)) -le $NREMAIN ]] ; then
    NREPEAT=$((${NREPEAT}+1))
fi
CNT=$(((${SGE_TASK_ID}-1)/$#))
BASEX=16
THISX=$((${BASEX}*${NREPEAT}))
THISPHS=$((${CNT}*${BASEX}))
#---------------------------------------------------------------------------
CAL=psa112_v005
CH=50_199_2
CAT=helm,misc,culgoora
#BIGSRCS=cyg,cas,vir,crab,her,hyd
#LITSRCS=`srclist.py -s 10/.150 -x 500/.150 --dec=20_50 --divstr=,`
#SRCS=$BIGSRCS,$LITSRCS
#SRCS=pic,hyd,for,J2214-170,cyg,crab,J1615-605,J1935-461,J2154-692,J2358-605
#SRCS=Sun,pic,hyd,for,J2214-170,cyg,crab,J1615-605,J1935-461,J2154-692,J2358-605,cen
SRCS=Sun,pic,hyd,for,J2214-170,cyg,crab,J1615-605,1932-464,J2154-692,J2358-605,cen
#ANTS=cross
#ANTS="(0,1,3,5,7,12,18,19,20,21,22,27,28,29)_(0,1,3,5,7,12,18,19,20,21,22,27,28,29),cross,-27_28,-0_12,-3_22"
ANTS="(1,5,19,20)_(0,1,3,5,7,12,18,19,20,21,22,27,28,29),cross,-27_28,-0_12,-3_22"
#ANTS="(1,5,7,18,19,20,21,28,29)_(1,5,7,18,19,20,21,28,29),cross"
#0,3,12,22,27
#ANTS="(22)_(0,12,27,1,5,7,18,19,20,21,28,29)"
#ANTS="(0,3,28,29,7,12,21,22)_(7,12,21,22),cross"
#ANTS="(0,3,28,29)_(0,3,28,29),cross"
#ANTS="-0_12,-3_22,-27_28"
#ANTS="(0,3,27)_(29,19,1,18,20),(0,3,27)_(0,3,27),cross"
#FITSRCS=`python -c "print '${SRCS}'.replace(',','/')"`
#PRMS="(${FITSRCS})=(jys)"
# good:.715
# all:.790, 0:.836, 1:.709, 3:.800, 5:.751, 7:.783, 12:.895, 18:.724, 19:.692
# 20:.733, 21:.780, 22:.836, 27:.859, 28:.778, 29:.728
# all-3bad:.749
#PRMS="(0/1/3/5/7/12/18/19/20/21/22/27/28/29)=(bp_r)"
PRMS="(1/5/19/20)=(phsoff_run6)"
#PRMS="(Sun/cen)=(jys/index/a1/a2/th)"
#PRMS="(Sun/pic/hyd/for/J2214-170/cyg/crab/J1615-605/1932-464/J2154-692/J2358-605/cen)=(jys/index)"
#PRMS="(27/18)=(amp)"
POL=yy

#---------------------------------------------------------------------------
echo "fitmdl.py -p $POL -C $CAL -c $CH -x $THISX --dphs=$THISPHS -s $SRCS -a \"$ANTS\" -P \"$PRMS\" $ARGS --baseport=54600 --cat=$CAT --master=\`qstat | qstat_to_hostport.py ${JOB_ID}\` | tee out.txt"
fitmdl.py -p $POL -C $CAL -c $CH -x $THISX --dphs=$THISPHS -s $SRCS -a $ANTS -P $PRMS --baseport=54600 --daemon=$SGE_TASK_ID $ARGS --cat=$CAT
