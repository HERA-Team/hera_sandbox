#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N img_zen
#$ -o grid_output/
#$ -l h_vmem=1G

SKIP=0
ENDCNT=200
CNT=0
FILE=`pull_args.py $*`
SNAP=168
#echo "mk_img.py --use_head --altmin=60 --skip=${SKIP} --cnt=${CNT} --facet_dec=-90_10  -c 144_370,475_900 --size=400 --res=0.4 -pxx -C psa331_v009_gc --fmt=${FILE}_%04d $FILE"
#mk_img.py --use_head --altmin=60 --skip=${SKIP} --cnt=${CNT} --facet_dec=-90_10  -c 144_370,475_900 --size=400 --res=0.4 -pxx -C psa331_v009_gc --fmt=${FILE}_%04d $FILE

echo "~/bin/mk_img.py -a cross,-24 --minuv=20 --altmin=30 --snap=${SNAP} -s zen --size=400 --res=0.4 -pxx -C psa455_v003_gc --fmt=${FILE}_%04d $FILE"
~/bin/mk_img.py -a cross,-24 --minuv=20 --altmin=30 -s zen --snap=${SNAP} --size=400 --res=0.4 -pxx -C psa455_v003_gc --fmt=${FILE}_%04d $FILE