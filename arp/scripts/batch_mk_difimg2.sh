#$ -S /bin/bash
POL=yy
CAL=pgb966_v003
ANTS=cross,-0_5,-0_6,-1_8,-3_11,-3_12,-4_15,-5_6,-12_13
#ANTS=cross
SZ=400
RES=.4
ALTMIN=45
#ALTMIN=30
UNIF=.01
NSIDE=512
C1=15
C2=185
dC=5
CHS=`python -c "for a in range($C1,$C2,$dC): print '%d_%d' % (a,a+5)"`
MYCHS=`pull_args.py $CHS`
for ch in $MYCHS ; do
    echo Working on channels: $ch
    FMT_FILE=pgb966_c${ch}_
    FMT=${FMT_FILE}%04d
    mk_difimg.py -p $POL -a $ANTS -C $CAL -c $ch --fmt=${FMT} --size=$SZ --res=$RES -o dim,dbm --altmin=$ALTMIN $* -u $UNIF
    cl_img.py -d cln --maxiter=10000 ${FMT_FILE}*.d[ib]m.fits
    mk_map.py --nside=$NSIDE -m ${FMT_FILE}map.fits ${FMT_FILE}*.bim.fits
done

