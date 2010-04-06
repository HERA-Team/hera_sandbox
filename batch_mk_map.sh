#$ -S /bin/bash
NSIDE=512
C1=15
C2=185
dC=5
CHS=`python -c "for a in range($C1,$C2,$dC): print '%d_%d' % (a,a+5)"`
MYCHS=`pull_args.py $CHS`
for ch in $MYCHS ; do
    echo Working on channels: $ch
    FMT_FILE=pgb966_c${ch}_
    #mk_map.py --nside=$NSIDE -m ${FMT_FILE}map.fits ${FMT_FILE}*[ab].bim.fits
    mk_map.py --nside=$NSIDE -m ${FMT_FILE}map.fits ${FMT_FILE}*.bim.fits
done
