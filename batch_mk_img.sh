#$ -S /bin/bash
POL=yy
CAL=pgb966_v002
SZ=400
RES=.4
ALTMIN=45
UNIF=.01
C1=15
C2=185
dC=5
CHS=`python -c "for a in range($C1,$C2,$dC): print '%d_%d' % (a,a+5)"`
MYCHS=`pull_args.py $CHS`
for ch in $MYCHS ; do
    echo Working on channels: $ch
    FMT_FILE=pgb966_c${ch}_
    FMT=${FMT_FILE}%04d
    echo '1/2 ---------------------------------------------------------'
    mk_img.py -p $POL -C $CAL -c $ch --fmt=${FMT}a --size=$SZ --res=$RES -o dim,dbm --altmin=$ALTMIN $* -x 2 --dphs=0 -u $UNIF
    echo '2/2 ---------------------------------------------------------'
    mk_img.py -p $POL -C $CAL -c $ch --fmt=${FMT}b --size=$SZ --res=$RES -o dim,dbm --altmin=$ALTMIN $* -x 2 --dphs=1 -u $UNIF
done

