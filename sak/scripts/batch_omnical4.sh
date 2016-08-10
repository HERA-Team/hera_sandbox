!!OUTDATED!!


#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N OMNICAL4
#$ -o grid_output/
#$ -l h_vmem=15G
#$ -t 1:20

echo 'here we go!'
source /usr/global/paper/CanopyVirtualEnvs/PAPER_Omni/bin/activate

FILES=`/data4/paper/2012EoR/pol_live/pull_args.py $*`

for FILE in $FILES; do
	arr=($(echo ${FILE} | tr "/" "\n"))
	echo omnical4.py ${FILE} -C psa6240_v003 -i redundantinfo_first_cal_2015_05_24_20_25_36.bin -r calpar_first_cal_2015_05_24_20_25_36.pmodel -o ${arr[0]} -u -s
	omnical4.py ${FILE} -C psa6240_v003 -i redundantinfo_first_cal_2015_05_24_20_25_36.bin -r calpar_first_cal_2015_05_24_20_25_36.pmodel -o ${arr[0]} -u -s
done
