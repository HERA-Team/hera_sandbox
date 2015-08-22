#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l paper
#$ -l h_vmem=2G
#$ -j y
#$ -N spectrum
#$ -o grid_output
shopt -s extglob

#SRCS=`cat ~/psa64_analysis/picA_spectrum/psa64_-45_dec_stripe.txt`
#SRCS=`cat ~/psa64_analysis/beamforms/arp_method/psa64_-12_dec_stripe.txt`
CALSRC=2331-416
#CALSRC=0213-132
SRCS=`cat ~/psa64_analysis/beamforms/psa64_pic_stripe_final.txt`
SRCS=`~/scripts/pull_args.py ${SRCS}`
for SRC in $SRCS
do
echo ~/scripts/beam_src_vs_ha.py \
 -s ${CALSRC} -C psa746_v012 -pxx --cat=southern_sky_v3 --sep_min=5 \
 lst*bm_${SRC}.F lst*bm_${CALSRC}.F
~/scripts/beam_src_vs_ha.py \
  -s ${CALSRC} -C psa746_v012 -pxx --cat=southern_sky_v3 --sep_min=5 \
 lst*bm_${SRC}.F  lst*bm_${CALSRC}.F
done
