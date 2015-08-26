#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N PSpecPrepU
#$ -o grid_output/
#$ -l h_vmem=5G
#$ -t 1:100

echo 'here we go!'
source activate PAPER

FILES=`python pull_args.py $*`

for FILE in $FILES; do
	echo python pspec_prep.py -C psa6240_v003 -p U -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${FILE}
	python pspec_prep.py -C psa6240_v003 -p U -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${FILE}
done
