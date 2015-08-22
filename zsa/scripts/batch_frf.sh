#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l paper
FILES=`pull_args.py $*`

for f in $FILES; do 
    echo fringe_rate_filter.py -a cross -p I -C psa6240_v003 --clean=1e-2 --max_fr_frac=1 --min_fr_frac=0 ${f}
    fringe_rate_filter.py -a cross -p I -C psa6240_v003 --clean=1e-2 --max_fr_frac=1 --min_fr_frac=0 ${f}
done

