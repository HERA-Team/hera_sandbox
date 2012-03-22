#$ -S /bin/bash
#---------------------------------------------------------------------------
# Admittedly complicated header to first split files between scripts, and if
# more nodes remain, split data within a file across nodes using decimation
# parameters.
ARGS=`pull_args.py -w $*`
echo $ARGS
NREPEAT=$(((${SGE_TASK_LAST}-${SGE_TASK_FIRST}+1)/$#))
NREMAIN=$(((${SGE_TASK_LAST}-${SGE_TASK_FIRST}+1)%$#))
if [[ $((${SGE_TASK_ID}%$#)) -le $NREMAIN ]] ; then
    NREPEAT=$((${NREPEAT}+1))
fi
CNT=$(((${SGE_TASK_ID}-1)/$#))
BASEX=8
THISX=$((${BASEX}*${NREPEAT}))
THISPHS=$((${CNT}*${BASEX}))
#---------------------------------------------------------------------------
CAL=psa112_v002_2
#CAL=pgb015_v004
#CAL=pgb050_v001
#_gps_v005
CH=60_200
#CH=200_250_5
#BIGSRCS=cyg,cas,vir,crab,her,hyd
BIGSRCS=pic,for
#LITSRCS=`srclist.py -s 10/.150 -x 500/.150 --dec=20_50 --ra=11_3 --divstr=,`
#SRCS=$BIGSRCS,$LITSRCS
SRCS=$BIGSRCS
POL=yy
#,Sun
#ANTS="(2,3,5,6)_(2,3,5,6)"
ANTS=cross,-1
#PRMS="(6)=(x/y/z/dly/off)"
#SPRMS=bp_r`
#PRMS="(0y/2y/3y/4y)=(amp)/1"
#SHPRMS="(2/3/4/5)=(amp)"
#PRMS="(2/3/4)=(amp)"
#PRMS="(2/3/4/5/6/8/9/10/11/12/14)=(dly/off)"
PRMS="(0/1/5/19/20/28/29)=(dly/off)"
#PRMS="(2/3/4)=dly/(12.85/1.65/9.81),(2/3/4)=off/0"
#7,13,15 outriggers, 1 is gom,0 is phase ref
#FITSRCS=`python -c "print '${SRCS}'.replace(',','/')"`
#PRMS="(${FITSRCS})=(jys)"
#PRMS="(${FITSRCS})=(jys)"
#PRMS="(Sun)=(jys/index/a1/a2)"
#---------------------------------------------------------------------------
#$ -o /home/jacobsda/jacobsda/gimage/psa112/grid_output
#$ -j y
#---------------------------------------------------------------------------
if ! ls /scratch/paper &> /dev/null; then
    mkdir /scratch/paper
fi
LOCAL=/scratch/paper/psa113/
if ! ls $LOCAL &> /dev/null; then
    echo Creating local data dir...
    mkdir $LOCAL
#    mv /tmp/lst*uv $LOCAL
fi
for F in $ARGS; do
    if ! ls ${LOCAL}${F} &> /dev/null;then
       echo preloading ${F}...  
       rsync -avz $SGE_O_HOST:${SGE_O_WORKDIR}/${F} ${LOCAL}
    fi
    if [ $? -ne 0 ]; then echo File transfer error!; exit; fi
    echo Finished preloading.
   TARGS=${LOCAL}`python -c "print '${F}'.split('/')[-1]"` ${TARGS}
   echo $TARGS
done
if  [ -z "${SHPRM+xxx}" ]; then
    echo Doing individual parameters
    echo "fitmdl.py -p $POL -C $CAL -c $CH -x $THISX --dphs=$THISPHS -s $SRCS -a \"$ANTS\" -P \"$PRMS\" $ARGS --master=\`qstat | qstat_to_hostport.py $JOB_ID\` | tee out.txt"
    fitmdl.py -p $POL -C $CAL -c $CH -x $THISX --dphs=$THISPHS -s $SRCS -a $ANTS -P $PRMS --daemon=$SGE_TASK_ID $TARGS
else
    echo Doing shard parameters
    echo "fitmdl.py -p $POL -C $CAL -c $CH -x $THISX --dphs=$THISPHS -s $SRCS -a \"$ANTS\" -S \"${SHPRMS}\" $ARGS  --master=\`qstat | qstat_to_hostport.py $JOB_ID\` | tee out.txt"
    fitmdl.py -p $POL -S ${SHPRMS} -C $CAL -c $CH -x $THISX --dphs=$THISPHS -s $SRCS -a $ANTS --daemon=$SGE_TASK_ID $TARGS
fi
