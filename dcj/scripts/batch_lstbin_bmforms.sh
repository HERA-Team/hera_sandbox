#$ -V
#$ -cwd
#$ -l paper
#$ -l h_vmem=2G
#$ -j y
#$ -N lstbin
#$ -o grid_output
shopt -s extglob

#SRCS=`cat ~/psa64_analysis/picA_spectrum/psa64_-45_dec_stripe.txt`
#SRCS=`cat ~/psa64_analysis/beamforms/arp_method/psa64_-12_dec_stripe.txt`
SRCS=`cat ~/psa64_analysis/beamforms/psa64_pic_stripe_final.txt`
SRCS=`~/scripts/pull_args.py ${SRCS}`
for SRC in $SRCS
do
echo ~/arp_scripts/lstbin.py --nogaps -C psa746_v011 --lst_res=85.9 ~/psa64_analysis/beamforms/July/*bm_${SRC} ~/psa64_analysis/beamforms/October/*bm_${SRC}
~/scripts/lstbin.py --nogaps -C psa746_v012 --lst_res=85.9 ~/psa64_analysis/beamforms/July/*bm_${SRC} ~/psa64_analysis/beamforms/October/*bm_${SRC}
echo ~/arp_scripts/fringe_rate_filter_static.py --maxfr=0.0001 --minfr=-0.0001 --clean=1e-6 lst*bm_${SRC}
~/scripts/fringe_rate_filter_static.py --maxfr=0.0001 --minfr=-0.0001 --clean=1e-6 lst*bm_${SRC}
done
