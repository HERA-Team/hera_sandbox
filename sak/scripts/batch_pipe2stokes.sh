#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N pipe2stokes
#$ -o grid_output/
#$ -l h_vmem=10G
#$ -t 1:100

echo 'here we go!'
source activate PAPER

FILES=`pull_args.py $*`

for FILE in $FILES; do
	echo python pipe2stokes.py ${FILE}
	python pipe2stokes.py ${FILE}
done
