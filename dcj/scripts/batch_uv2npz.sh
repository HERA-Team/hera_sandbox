#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l paper
#$ -l h_vmem=2G
#$ -j y
#$ -N npz
#$ -o grid_output
CAL=psa5842_v001
#SRCS=`cat psa64_-45_dec_stripe.txt`
SRCS=`cat ~/psa64_analysis/beamforms/psa64_pic_stripe_final.txt`
SRCs=`~/scripts/pull_args.py ${SRCS}`
for SRC in $SRCs
do
echo ~/scripts/beamform_uv2npz.py -C $CAL -p xx  -s $SRC --cat=southern_sky_v3 \
~/psa64_analysis/beamforms/LSTbinned/*bm_${SRC}.F

~/scripts/beamform_uv2npz.py -C $CAL -p xx  -s $SRC --cat=southern_sky_v3 \
 ~/psa64_analysis/beamforms/LSTbinned/*bm_${SRC}.F

done

