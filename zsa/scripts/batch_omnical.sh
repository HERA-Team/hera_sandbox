#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -e grid_output
#$ -o grid_output
#$ -l paper
#$ -l h_vmem=16G
DIRS=`pull_args.py $*`
CALFILE=psa6240_v003
EXT=uvcRREcAC

for dir in ${DIRS}; do
    TAG1=`cut -d /  -f 6`
    TAG2=`cut -d /  -f 7 | cut -d . -f 2`
    for CHUNK in `seq .1 .1 .6`; do
        echo omnical_PSA64.py -C ${CALFILE} -p xx,yy --tag=${TAG1}_${TAG2}.${CHUNK} $dir/zen.*.${CHUNK}*.${EXT}
        omnical_PSA64.py -C ${CALFILE} -p xx,yy --tag=${TAG1}_${TAG2}.${CHUNK} $dir/zen.*.${CHUNK}*.${EXT}
    done;
done;
    
