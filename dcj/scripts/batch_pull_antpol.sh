#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=2G
#$ -j y
#$ -N pull
#$ -o grid_output
#$ -q all.q
CAL=psa6622_v002
SEPS=0,2
POL=xx
PATH=/home/jacobsda/src/anaconda/bin/:$PATH
source activate PAPER
FILES=`~/scripts/pull_args.py $*`
echo pulling data from $FILES
OUTRIGGERS="64_112,64_113,64_114,64_115,64_116,64_117,64_118,64_119,64_120,64_121,64_122,64_123,64_124,64_125,64_126,64_127,49_112,49_113,49_114,49_115,49_116,49_117,49_118,49_119,49_120,49_121,49_122,49_123,49_124,49_125,49_126,49_127"

for SEP in $SEPS
do
mkdir sep$SEP
cd sep$SEP
BLS=$(grid2ant.py -C ${CAL} --seps=${SEP})
~/scripts/pull_antpols.py -a ${BLS} $FILES
cd ..
done
