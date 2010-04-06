#$ -S /bin/bash
POL=yy
#CAL=pgb015_v005
CAL=pgb564_ver005
#ANTS=cross,-0_5,-0_6,-5_6,-3_12
ANTS=cross
SZ=400
RES=.4
ALTMIN=45
CEN="11:20_40:00"
NSIDE=512
C1=60
C2=180
dC=5
CHS=`python -c "for a in range($C1,$C2,$dC): print '%d_%d' % (a,a+5)"`
MYCHS=`pull_args.py $CHS`
for ch in $MYCHS ; do
    echo Working on channels: $ch
    FMT_FILE=${CAL}_c${ch}_
    FMT=${FMT_FILE}%04d
    mk_img.py -s $CEN -p $POL -a $ANTS -C $CAL -c $ch --fmt=${FMT}a --size=$SZ --res=$RES -o dim,dbm --altmin=$ALTMIN -x 2 --dphs=0 $*
    mk_img.py -s $CEN -p $POL -a $ANTS -C $CAL -c $ch --fmt=${FMT}b --size=$SZ --res=$RES -o dim,dbm --altmin=$ALTMIN -x 2 --dphs=1 $*
    cl_img.py -d cln --maxiter=500 --div --tol=1e-7 -r radial -o bim,rim ${FMT_FILE}*.d[ib]m.fits
done

