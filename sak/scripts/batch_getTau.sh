#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N TAUXY
#$ -o grid_output/
#$ -l h_vmem=5G
#$ -t 1:20

echo 'here we go!'
source /usr/global/paper/CanopyVirtualEnvs/PAPER_Omni/bin/activate

#call as psa????
FOLDERS=`/data4/paper/2012EoR/pol_live/pull_args.py $*`

for FOLDER in $FOLDERS; do
	echo python getTau.py ${FOLDER} 0 redundantinfo_first_cal_2015_06_11_12_29_37.bin
	python getTau.py ${FOLDER} 0 redundantinfo_first_cal_2015_06_11_12_29_37.bin
	echo python getTau.py ${FOLDER} 15 redundantinfo_first_cal_2015_06_11_12_29_37.bin
	python getTau.py ${FOLDER} 15 redundantinfo_first_cal_2015_06_11_12_29_37.bin
	echo python getTau.py ${FOLDER} 30 redundantinfo_first_cal_2015_06_11_12_29_37.bin
	python getTau.py ${FOLDER} 30 redundantinfo_first_cal_2015_06_11_12_29_37.bin
	echo python getTau.py ${FOLDER} 45 redundantinfo_first_cal_2015_06_11_12_29_37.bin
	python getTau.py ${FOLDER} 45 redundantinfo_first_cal_2015_06_11_12_29_37.bin
done
