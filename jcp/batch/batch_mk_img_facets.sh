#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N img_facets
#$ -o grid_output/

ANTS="cross,-16,-31,-3_13,-3_14,-3_29,-4_14,-4_15,-4_30,-4_32,-5_15,-5_32,-6_32,-13_20,-14_21,-15_21,-15_22,-15_23,-20_29,-21_30,-21_32,-22_32,-23_32"

SKIP=0
ENDCNT=200
CNT=0
FILE=`pull_args.py $*`

echo "dj_mk_img.py --use_head -a $ANTS --altmin=30 --skip=${SKIP} --cnt=${CNT} --facet_dec=-10_90  -c 80_180 --size=400 --res=0.4 -pxx -C pgb322_v008_gc --fmt=${FILE}_%s_%04d $FILE"
dj_mk_img.py --use_head -a $ANTS --altmin=30 --cnt=${CNT} --facet_dec=-10_90  -c 80_180 --size=400 --res=0.4 -pxx -C pgb322_v008_gc --fmt=${FILE}_%s_%04d $FILE