#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N FINE_IMG
#$ -o grid_output/
#$ -l h_vmem=1G

CAL=psa819_v007
FILES=`pull_args.py $*`

FSTART=150
FSTOP=450
NCHAN=$(($FSTOP-$FSTART))
df=$(($NCHAN/32))

for FILE in $FILES; do
    echo ===============================
    for i in `seq 0 31`; do
        echo =========
        f0=$(($i*$df+$FSTART))
        f1=$(($f0+$df))
    
        for pol in xx xy yx yy; do
            echo mk_img.py -C ${CAL} -s zen --size=500 --res=0.3 -c ${f0}_${f1} -a cross -p ${pol} ${FILE} \
                --fmt=${FILE}_FB${i}${pol}_%04d
            mk_img.py -C ${CAL} -s zen --size=500 --res=1. -c ${f0}_${f1} -a cross -p ${pol} ${FILE} \
                --fmt=${FILE}_FB${i}${pol}_%04d
        
        done    
    done
done
