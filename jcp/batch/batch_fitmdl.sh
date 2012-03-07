#$ -S /bin/bash
#---------------------------------------------------------------------------
# Admittedly complicated header to first split files between scripts, and if
# more nodes remain, split data within a file across nodes using decimation
# parameters.
source .bashrc
#echo $PATH

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
CAL=pgb322_v003_gc
CH=160_170
#CAT=helm,misc,culgoora

SRCS=cyg,cas
ANTS=cross,-16,-3_14
PRMS="(0/2/3/4/5/6/8/9/10/11/12/13/14/15/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32)=(x/y/z/off)"
POL=xx

#---------------------------------------------------------------------------
echo "fitmdl.py -p $POL -C $CAL -c $CH -s $SRCS -a \"$ANTS\" -P \"$PRMS\" -x $THISX --dphs=THISPHS $ARGS --baseport=53000 --master=`qstat | qstat_to_hostport.py ${JOB_ID}` | tee /data1/paper/pober/out.txt"
#echo "master=`qstat | qstat_to_hostport.py ${JOB_ID}`"
fitmdl.py -p $POL -C $CAL -c $CH -s $SRCS -a $ANTS -P $PRMS -x $THISX --dphs=$THISPHS --baseport=53000 --daemon=$SGE_TASK_ID $ARGS
