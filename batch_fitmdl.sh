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
BASEX=2
THISX=$((${BASEX}*${NREPEAT}))
THISPHS=$((${CNT}*${BASEX}))
#---------------------------------------------------------------------------
CAL=psa112_v003
CH=40_205
#BIGSRCS=cyg,cas,vir,crab,her,hyd
#LITSRCS=`srclist.py -s 10/.150 -x 500/.150 --dec=20_50 --divstr=,`
#SRCS=$BIGSRCS,$LITSRCS
SRCS=pic,hyd,for,J2214-170,cyg,crab,J1615-605,J1935-461,J2154-692,J2358-605
#SRCS=Sun
ANTS=cross
#FITSRCS=`python -c "print '${SRCS}'.replace(',','/')"`
#PRMS="(${FITSRCS})=(jys)"
PRMS="(J2358-605)=(jys/index/ra/dec)"
POL=yy

#---------------------------------------------------------------------------
echo "fitmdl.py -p $POL -C $CAL -c $CH -x $THISX --dphs=$THISPHS -s $SRCS -a \"$ANTS\" -P \"$PRMS\" $ARGS --baseport=54600 --master=\`qstat | qstat_to_hostport.py ${JOB_ID}\` | tee out.txt"
fitmdl.py -p $POL -C $CAL -c $CH -x $THISX --dphs=$THISPHS -s $SRCS -a $ANTS -P $PRMS --baseport=54600 --daemon=$SGE_TASK_ID $ARGS
