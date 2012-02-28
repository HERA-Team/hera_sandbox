#$ -S /bin/bash
POL=yy
CAL=pgb966_v003
SZ=400
RES=.4
ALTMIN=45
#C1=15
#C2=185
#dC=5
C1=150
C2=155
dC=1
SRC="12:00_40:00"
CHS=`python -c "for a in range($C1,$C2,$dC): print '%d_%d' % (a,a+${dC})"`
MYCHS=`pull_args.py $CHS`
for ch in $MYCHS ; do
    echo Working on channels: $ch
    FMT_FILE=pspec1_pgb966_c${ch}_
    FMT=${FMT_FILE}%04d
    echo '1/2 ---------------------------------------------------------'
    mk_img.py -s $SRC -p $POL -C $CAL -c $ch --fmt=${FMT}a --size=$SZ --res=$RES -o dim,dbm --altmin=$ALTMIN $* -x 2 --dphs=0
    echo '2/2 ---------------------------------------------------------'
    mk_img.py -s $SRC -p $POL -C $CAL -c $ch --fmt=${FMT}b --size=$SZ --res=$RES -o dim,dbm --altmin=$ALTMIN $* -x 2 --dphs=1
done

