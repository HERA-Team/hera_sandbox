#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N img_sngl
#$ -o img_sngl


FILES=`pull_args.py $*`

echo "mk_img.py -s cen --altmin=60  -c 400_600 --size=400 --res=0.4 -pxx -C psa331_v008_gc $FILES --fmt=${FILES}_${SGE_TASK_ID}_%02f $FILES"
mk_img.py -s cen --altmin=60  -c 400_600 --size=400 --res=0.4 -pxx -C psa331_v008_gc $FILES --fmt=${FILES}_${SGE_TASK_ID}_%02f $FILES
