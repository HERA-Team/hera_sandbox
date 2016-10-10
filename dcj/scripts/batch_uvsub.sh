#$ -S /bin/bash
#$ -V
#$ -j y
#$ -cwd
#$ -o grid_output/
#$ -N filter_sun
CAL=psa331_v008_gc


ARGS=`pull_args.py $*`
for FILE in ${ARGS}
do 
    echo uv_addsub.py --sub $FILE ${FILE}.e_Sun
    uv_addsub.py --sub $FILE ${FILE}.e_Sun
done
#echo filter_src.py -e -C $CAL -s $SRCS -d 5 -r 5 --clean=5e-4 $ARGS
#filter_src.py -e -C $CAL -s $SRCS -d 5 -r 5 --clean=5e-4 $ARGS