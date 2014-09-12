# !/bin/bash

#PBS -S /bin/bash
#PBS -N monte_carlo_analysis_3_carver
#PBS -j eo
#PBS -l nodes=1:ppn=1,walltime=00:05:00,pvmem=20GB
#PBS -q regular
#PBS -A m1871


codeLoc="/global/homes/a/acliu/globalSig/fq_120_150_testCase"

module load python/2.7.3
module load numpy
module load matplotlib

echo "******************************************************************"
echo "Here's the script that was run"
echo `cat $0`
echo "******************************************************************"

cd $codeLoc

echo "Analyzing MC results"
date
python "$codeLoc/monte_carlo_analysis_3.py"
python "$codeLoc/create_analysis_directories.py"
date
echo "...all done!"