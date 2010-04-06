#$ -S /bin/bash
POL=yy
#CAL=pgb015_v005
CAL=pgb564_ver006
ANTS=cross,-0_5,-0_6,-5_6,-3_12
#ANTS=cross
#ANTS=cross,-1
SZ=400
RES=.4
ALTMIN=45
NSIDE=512
C1=60
C2=180
dC=5
CHS=`python -c "for a in range($C1,$C2,$dC): print '%d_%d' % (a,a+5)"`
MYCHS=`pull_args.py $CHS`
for ch in $MYCHS ; do
  for cen in "11:20_40:00" "11:20_30:00" "12:00_40:00" "12:00_30:00"; do
  #for cen in "12:00_40:00"; do
    echo Working on channels: $ch, $cen
    CEN_FMT=`python -c "print '${cen}'.replace(':','')"`
    FMT_FILE=${CEN_FMT}_c${ch}_
    FMT=${FMT_FILE}%02d
    mk_img.py -s $cen -p $POL -a $ANTS -C $CAL -c $ch --fmt=${FMT}a --size=$SZ --res=$RES -o dim,dbm --altmin=$ALTMIN -x 2 --dphs=0 $*
    mk_img.py -s $cen -p $POL -a $ANTS -C $CAL -c $ch --fmt=${FMT}b --size=$SZ --res=$RES -o dim,dbm --altmin=$ALTMIN -x 2 --dphs=1 $*
    #cl_img.py -d cln --maxiter=500 --div --tol=1e-7 -r radial -o bim,rim ${FMT_FILE}*.d[ib]m.fits
  done
done

