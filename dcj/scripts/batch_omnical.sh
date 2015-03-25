#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l paper
#$ -l h_vmem=8G
#$ -j y
#$ -N omnical
#$ -o grid_output
CAL=psa6942_v000
#NOTE: To be run in ~jacobsda/storage/psa128
canopy-PAPER_Omni
FILES=`~/scripts/pull_args.py $*`
for FILE in $FILES
do
echo Processing $FILE
echo omnical_PSA128.py -C ${CAL} -pxx -o 2014_epoch3/ -i omnicals/redundantinfo_first_cal_psa128_6942.bin -r omnicals/calpar_first_cal_psa128_6942.p -u -d psa128_6942 -t O -s $FILE --add
time omnical_PSA128.py -C ${CAL} -pxx -o 2014_epoch3/ -i omnicals/redundantinfo_first_cal_psa128_6942.bin -r omnicals/calpar_first_cal_psa128_6942.p -u -d psa128_6942 -t O -s $FILE --add

done

