#$ -S /bin/bash
POL=yy
CAL=psa112_v003
ANTS=cross
SZ=400
RES=.4
#ALTMIN=45
ALTMIN=30
NSIDE=512
C1=40
C2=205
dC=5
MINUV=20
CHS=`python -c "for a in range($C1,$C2,$dC): print '%d_%d' % (a,a+5)"`
MYCHS=`pull_args.py $CHS`
for ch in $MYCHS ; do
    echo Working on channels: $ch
    FMT_FILE=${CAL}_c${ch}_
    FMT=${FMT_FILE}%04d
    mk_img.py -p $POL -a $ANTS -C $CAL -c $ch --fmt=${FMT} --size=$SZ --res=$RES -o dim,dbm --altmin=$ALTMIN $*
    cl_img.py -d cln --maxiter=500 --div --tol=1e-7 -r radial -o bim,rim ${FMT_FILE}*.d[ib]m.fits
    mk_map.py --nside=$NSIDE -m ${FMT_FILE}bmap.fits ${FMT_FILE}*.bim.fits
    mk_map.py --nside=$NSIDE -m ${FMT_FILE}rmap.fits ${FMT_FILE}*.rim.fits
done

