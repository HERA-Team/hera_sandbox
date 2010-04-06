#$ -S /bin/bash
POL=yy
CAL=pgb564_ver006
ANTS=cross,-0_5,-0_6,-5_6,-3_12
SZ=400
RES=.4
ALTMIN=45
NSIDE=512
C1=60
C2=180
dC=5
MINUV=0
CHS=`python -c "for a in range($C1,$C2,$dC): print '%d_%d' % (a,a+5)"`
MYCHS=`pull_args.py $CHS`
for ch in $MYCHS ; do
    echo Working on channels: $ch
    FMT_FILE=${CAL}_c${ch}_
    FMT=${FMT_FILE}%04d
    mk_img.py -a $ANTS -p $POL -C $CAL -c $ch --fmt=${FMT}a --size=$SZ --res=$RES -o dim,dbm --altmin=$ALTMIN --minuv=$MINUV -x 2 --dphs=0 $*
    mk_img.py -a $ANTS -p $POL -C $CAL -c $ch --fmt=${FMT}b --size=$SZ --res=$RES -o dim,dbm --altmin=$ALTMIN --minuv=$MINUV -x 2 --dphs=1 $*
    cl_img.py -d cln --maxiter=500 --div --tol=1e-7 -r radial -o bim ${FMT_FILE}*.d[ib]m.fits
    mk_map.py --nside=$NSIDE -m ${FMT_FILE}bmap_a.fits ${FMT_FILE}*a.bim.fits
    mk_map.py --nside=$NSIDE -m ${FMT_FILE}bmap_b.fits ${FMT_FILE}*b.bim.fits
done

