#$ -S /bin/bash
#$ -V
#$ -j y
#$ -cwd
#$ -o grid_output/
#$ -l h_vmem=100M

FILES=`pull_args.py $*`
for FILE in $FILES
do
echo "cl_img.py -d cln -r radial ${FILE} --maxiter=2000 --tol=1e-4"
cl_img.py -d cln -r radial ${FILE} --maxiter=2000 --tol=1e-4
done
#C1=15
#C2=185
#dC=5
#CHS=`python -c "for a in range($C1,$C2,$dC): print '%d_%d' % (a,a+5)"`
#MYCHS=`pull_args.py $CHS`
#for ch in $MYCHS ; do
#    echo Working on channels: $ch
#    FMT_FILE=pgb966_c${ch}_
#    cl_img.py -d cln --maxiter=10000 ${FMT_FILE}*.d[ib]m.fits
#done
