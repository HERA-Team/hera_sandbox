#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N img_psa64
#$ -o grid_output/
#$ -e grid_output/
POL=yy
CAL=psa746_wacky
SZ=400
RES=.4
ALTMIN=30
#C1=0
#C2=512
C1=819
#C2=921
#C1=921
C2=1331
#C2=1536
dC=5

CHS=`python -c "for a in range($C1,$C2,$dC): print '%d_%d' % (a,a+${dC}-1)"`
MYCHS=`pull_args.py $CHS`
for ch in $MYCHS; do
    echo Working on channels: $ch
    FMT_FILE=psa747_wacky_c${ch}_
    FMT=${FMT_FILE}%04d
    mk_img.py -p $POL -a cross -C $CAL -c $ch --fmt=$FMT --size=$SZ --res=$RES -o dim,dbm --altmin=$ALTMIN -s "21:52:10.67_-30:43:17.5" --zen_phased --skip_amp $*
    #mk_img.py -p $POL -a cross -C $CAL -c $ch --fmt=$FMT --size=$SZ --res=$RES -o dim,dbm --altmin=$ALTMIN -s "12:00:00_-30:43:17.5" --zen_phased --skip_amp $*
done

