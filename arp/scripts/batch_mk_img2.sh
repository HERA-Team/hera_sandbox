#$ -S /bin/bash
POL=yy
CAL=pgb015_v005
#ANTS="(0,1,3,5,7,12,18,19,20,21,22,27,28,29)_(0,1,3,5,7,12,18,19,20,21,22,27,28,29),cross,-27_28,-0_12,-3_22"
ANTS="cross,-1,-7,-13,-15"
#ANTS=cross
SZ=400
RES=.4
ALTMIN=30
NSIDE=512
C1=240
C2=720
dC=20
MINUV=20
CHS=`python -c "for a in range($C1,$C2,$dC): print '%d_%d' % (a,a+$dC-1)"`
MYCHS=`pull_args.py $CHS`
for ch in $MYCHS ; do
    echo Working on channels: $ch
    FMT_FILE=${CAL}_c${ch}_
    FMT=${FMT_FILE}%04d
    mk_img.py -p $POL -a $ANTS -C $CAL -c $ch --fmt=${FMT} --size=$SZ --res=$RES -o dim,dbm --altmin=$ALTMIN --minuv=$MINUV $*
    cl_img.py -d cln --maxiter=500 --div --tol=1e-7 -r radial -o bim ${FMT_FILE}*.d[ib]m.fits
    mk_map.py --nside=$NSIDE -m ${FMT_FILE}bmap.fits ${FMT_FILE}*.bim.fits
    #mk_map.py --nside=$NSIDE -m ${FMT_FILE}rmap.fits ${FMT_FILE}*.rim.fits
done

