#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N BYEXTALK
#$ -o grid_output/
#$ -l h_vmem=5G
#$ -t 1:50

echo 'here we go!'
source /usr/global/paper/CanopyVirtualEnvs/PAPER_Omni/bin/activate
#call psa????
FOLDERS=`/data4/paper/2012EoR/pol_live/pull_args.py $*`

for FOLDER in $FOLDERS; do
	echo python xtalk5.py -C psa6240_v003 -s ${FOLDER} ${FOLDER}/*.uvcRREcC
	python xtalk5.py -C psa6240_v003 -s ${FOLDER} ${FOLDER}/*.uvcRREcC
done
