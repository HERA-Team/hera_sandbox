#! /usr/bin/env bash

CH_LIST="240 360 480 600 719"

for FILE in $*; do
    for CH in $CH_LIST; do
        TAG=`python -c "print '${FILE}'[len('zen.2455'):len('zen.2455')+6]"`
        TAG=${TAG}`python -c "print '${FILE}'[-3:]"`
        TAG=${TAG}_c${CH}
        echo $TAG
        mk_img.py -a "cross,-1,-7,-13,-15" -p yy -C pgb015_v005 -s 13:30_40:00 -c $CH --size=400 --res=.4 --cnt=0 ${FILE}
        mv im0000.dim.fits 1330_4000_${TAG}.dim.fits
        mv im0000.dbm.fits 1330_4000_${TAG}.dbm.fits
        cl_img.py -d cln --div -r natural --maxiter=500 1330_4000_${TAG}.d[ib]m.fits -o bim,rim --tol=1e-7
    done
done
