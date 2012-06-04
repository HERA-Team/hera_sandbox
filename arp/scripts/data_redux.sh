#! /bin/bash

xrfi_simple.py -c 0_130,755_777,1540,1704,1827,1868,1885_2047 --df=6 *.uvGC
# XXX is the "pol" argument below necessary?
ddr_filter_coarse.py -p xx --clean=1e-3 --maxbl=300 --invert *.uvGCR
xrfi_simple.py -n 4 *.uvGCRE
for FILE in `ls -d *.uvGCR`; do
    echo Copying ${FILE}ER/flags to ${FILE}
    cp ${FILE}/flags ${FILE}/flags_bk
    cp ${FILE}ER/flags ${FILE}/flags
    #rm -rf ${FILE}[DEF]
done
#ddr_filter_coarse.py -p xx --clean=1e-4 --maxbl=300 *.uvGCR
#rm -rf *.uvGCRER
#    #if ! ls ${FILE}drxam &> /dev/null; then
#    #fi
