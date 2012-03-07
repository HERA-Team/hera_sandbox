#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N img_snapshot
#$ -o grid_output/
POL=xx
CAL=pgb322_v008_gc
ANTS="cross,-16,-31,-3_13,-3_14,-3_29,-4_14,-4_15,-4_30,-4_32,-5_15,-5_32,-6_32,-13_20,-14_21,-15_21,-15_22,-15_23,-20_29,-21_30,-21_32,-22_32,-23_32"
SZ=400
RES=.4
ALTMIN=30
NSIDE=512
MINUV=20
SNAP=120

ARGS=`pull_args.py $*`
for FILE in $ARGS; do
    echo Working on file: $FILE
    FMT_FILE=${FILE}_snap%04d
    mk_img.py -p $POL -a $ANTS -C $CAL -c 80_180 --fmt=$FMT_FILE --size=$SZ --res=$RES -o dim,dbm --snap=$SNAP --altmin=$ALTMIN --minuv=$MINUV -s zen $FILE
    cl_img.py -d cln --maxiter=10000 --div --tol=1e-6 -r radial -o bim ${FILE}_snap*.d[ib]m.fits
done

