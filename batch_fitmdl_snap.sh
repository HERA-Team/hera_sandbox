#$ -S /bin/bash
#---------------------------------------------------------------------------
# For fitmdl in snap mode, need to split by ch
NCHAN=1024
CH=200_302_2
MYCH=`python -c "import aipy; print ' '.join(map(str, aipy.scripting.parse_chans('${CH}', $NCHAN)))"`
echo $MYCH
MYCH=`pull_args.py -w ${MYCH}`
echo $MYCH
MYCH=`python -c "print '${MYCH}'.replace(' ',',')"`
echo $MYCH
#---------------------------------------------------------------------------
ARGS=$*
PORT=54600
MAXITER=250
XTOL=1e-2
FTOL=1e-2
#---------------------------------------------------------------------------
CAL=pgb015_v005
CAT=helm,misc,culgoora
SRCS=cyg,cas,vir,crab,Sun
ANTS=cross,-1,-7,-13,-15
#FITSRCS=`python -c "print '${SRCS}'.replace(',','/')"`
FITSRCS=cyg/cas/vir/crab/Sun
PRMS="(${FITSRCS})=(jys)"
PRMS="${PRMS},(0/2/3/4/5/6/8/9/10/11/12/14)=(amp)"
POL=yy

#---------------------------------------------------------------------------
echo "fitmdl.py -p $POL -C $CAL -c $CH -s $SRCS -a $ANTS -P \"$PRMS\" --baseport=$PORT --cat=$CAT --master=\`qstat | qstat_to_hostport.py ${JOB_ID}\` --snap --remem --maxiter=$MAXITER --xtol=$XTOL --ftol=$FTOL -q $ARGS | tee out.txt"
echo "fitmdl.py -p $POL -C $CAL -c $MYCH -s $SRCS -a $ANTS -P \"$PRMS\" --baseport=$PORT --cat=$CAT --daemon=$SGE_TASK_ID --snap $ARGS"
fitmdl.py -p $POL -C $CAL -c $MYCH -s $SRCS -a $ANTS -P $PRMS --baseport=$PORT --cat=$CAT --daemon=$SGE_TASK_ID --snap $ARGS
