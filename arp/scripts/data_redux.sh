#! /bin/bash

UVR=$*

#xrfi_simple.py -c 0_130,755_777,1540,1704,1827,1868,1885_2047 --df=6 $UV
# XXX is the "pol" argument in ddr_filter_coarse necessary?
for FILE in $UVR; do
    TRIPLET=`get_uv_neighbor.py $FILE`
    if ! ls ${FILE}/flags_bk &> /dev/null; then
        ddr_filter_coarse.py -p xx -a all --clean=1e-3 --maxbl=300 --output=ddr --invert $TRIPLET
        xrfi_simple.py -n 4 ${FILE}E
        rm -rf ${FILE}E
        echo Copying ${FILE}ER/flags to ${FILE}
        cp ${FILE}/flags ${FILE}/flags_bk
        cp ${FILE}ER/flags ${FILE}/flags
        #rm -rf ${FILE}ER
    fi
    ddr_filter_coarse.py -p xx -a all --clean=1e-4 --maxbl=300 $TRIPLET
done
