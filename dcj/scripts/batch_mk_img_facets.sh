#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N img_facets
#$ -o grid_output/
#$ -l h_vmem=1G

SKIP=0
ENDCNT=200
CNT=0
FILE=`pull_args.py $*`

#echo "mk_img.py --use_head --altmin=60 --skip=${SKIP} --cnt=${CNT} --facet_dec=-90_10  -c 144_370,475_900 --size=400 --res=0.4 -pxx -C psa331_v009_gc --fmt=${FILE}_%04d $FILE"
#mk_img.py --use_head --altmin=60 --skip=${SKIP} --cnt=${CNT} --facet_dec=-90_10  -c 144_370,475_900 --size=400 --res=0.4 -pxx -C psa331_v009_gc --fmt=${FILE}_%04d $FILE

echo "mk_img.py -a cross,-24 --minuv=20 --use_head --altmin=30 --skip=${SKIP} --cnt=${CNT} --facet_dec=-90_10  --size=400 --res=0.4 -pxx -C psa455_v003_gc --use_facet_num --fmt=${FILE}_%s_%04d $FILE"
mk_img.py -a cross,-24 --minuv=20 --use_head --altmin=30 --skip=${SKIP} --cnt=${CNT} --facet_dec=-90_10 --size=400 --res=0.4 -pxx -C psa455_v003_gc --use_facet_num --fmt=${FILE}_%s_%04d $FILE