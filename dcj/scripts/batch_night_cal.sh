#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l paper
#$ -l h_vmem=8G
#$ -j y
#$ -N nightcal
#$ -o grid_output
shopt -s extglob
. /usr/global/paper/paperenv.sh
workon dcj-PAPER
CALNIGHT=2455962
NIGHTS=`~/scripts/file_date.py -u zen.!($CALNIGHT).?????.uv*T`
NIGHTS=`~/scripts/pull_args.py ${NIGHTS}`
for NIGHT in $NIGHTS
do
echo  ~/scripts/night_cal.py  zen.$CALNIGHT.*T zen.$NIGHT.*T -C psa898_v003  --bin=0.02   --antcal  --apply --plot
 ~/scripts/night_cal.py  zen.$CALNIGHT.*T zen.$NIGHT.*T -C psa898_v003  --bin=0.02   --antcal  --apply --plot
done
