#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=3G
#$ -N MK_STOKES
source activate PAPER
ARGS=`pull_args.py $*` #hand this script one polarization
PATH2CAPO='/home/saulkohn/ReposForCanopy/capo/' 
STOKES="I,Q,U,V"

for f in ${ARGS}; do
    fS=""
    wd=`cut -d "z" -f 1 <<< "$f"`
    fQ=z`cut -d "z" -f 2 <<< "$f"`
    p=${fQ:18:2}
    for pp in xx xy yx yy; do
        fS+="${f/$p/$pp} "
    done
    echo ${PATH2CAPO}/pspec_pipeline/mk_stokes.py --stokes=${STOKES} $fS
    ${PATH2CAPO}/pspec_pipeline/mk_stokes.py --stokes=${STOKES} $fS
done   
