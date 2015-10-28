#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N byebyeRFI
#$ -o grid_output/
#$ -l h_vmem=5G
#$ -t 1:100

echo 'here we go!'
source activate PAPER

FILES=`python pull_args.py $*`

for FILE in $FILES; do
	echo xrfi_simple.py -n 3 ${FILE} 
	xrfi_simple.py -n 3 ${FILE} 
done
