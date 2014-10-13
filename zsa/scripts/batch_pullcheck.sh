#$ -S /bin/bash
#$ -V
#$ -cwd

DIRS=`pull_args.py $*`
for DIR in $DIRS; do
    for file in $DIR/*.uvcRREc; do
        echo checking ${file}A
        python -c "import aipy as a; uv = a.miriad.UV('${file}A')" || (rm -rf ${file}A; pull_antpols.py -p xx,yy -a all $file)
    done
done
