#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N OMNI4_EP2
#$ -o grid_output/
#$ -l h_vmem=15G
#$ -t 1:80

echo 'This is omnical4 for PSA-64 Epoch 2'
echo 'Valid for 2456310 onwards [psa63{1,2,3,4,5,6,7}?]'
echo 'here we go!'
source /usr/global/paper/CanopyVirtualEnvs/PAPER_Omni/bin/activate

FILES=`/data4/paper/2012EoR/pol_live/pull_args.py $*`

for FILE in $FILES; do
	arr=($(echo ${FILE} | tr "/" "\n"))
	echo omnical4.py ${FILE} -C psa6240_v003 -i redundantinfo_first_cal_2015_06_11_13_27_11.bin -r calpar_first_cal_2015_06_11_13_27_11.pmodel -o ${arr[0]} -u -s --skip_sun
	omnical4.py ${FILE} -C psa6240_v003 -i redundantinfo_first_cal_2015_06_11_13_27_11.bin -r calpar_first_cal_2015_06_11_13_27_11.pmodel -o ${arr[0]} -u -s --skip_sun
done
