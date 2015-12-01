#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N APPLYCAL
#$ -o grid_output/
#$ -l h_vmem=2G
#$ -t 1:20

echo 'here we go!'
source /usr/global/paper/CanopyVirtualEnvs/PAPER_Omni/bin/activate

FILES=`/data4/paper/2012EoR/pol_live/pull_args.py $*`

for FILE in $FILES; do
	echo apply_cal.py ${FILE} -C psa6240_v003
	apply_cal.py ${FILE} -C psa6240_v003
done
