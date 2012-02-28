#$ -S /bin/bash
POL=yy
CAL=pgb564_ver006
CENS="11:20_30:00 11:20_40:00 12:00_30:00 12:00_40:00"
ANTS=cross
SZ=400
RES=.4
ALTMIN=60
CHS=`python -c "print ' '.join(['%d_%d' % (c,c+4) for c in range(60,180,5)])"`
MYCHS=`pull_args.py $CHS`

for ch in $MYCHS; do
    echo Working on channels: $ch
    for cen in $CENS; do
        echo Working on center: $cen
        FMT_FILE=`python -c "print 'tsys_${cen}_c${ch}_'.replace(':','')"`
        FMT=${FMT_FILE}%02d
        #mk_img.py -p $POL -C $CAL -a $ANTS -c $ch -s $cen --fmt=${FMT}a --size=$SZ --res=$RES -o dim,dbm --altmin=$ALTMIN $* -x 2 --dphs=0 $*
        #mk_img.py -p $POL -C $CAL -a $ANTS -c $ch -s $cen --fmt=${FMT}b --size=$SZ --res=$RES -o dim,dbm --altmin=$ALTMIN $* -x 2 --dphs=1 $*
        cl_img.py -d cln --maxiter=500 --div --tol=1e-7 -o bim,rim ${FMT_FILE}*.d[ib]m.fits
    done
done
