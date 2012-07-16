#! /bin/bash

UV=$*
UV_R=`python -c "print ' '.join(map(lambda x: x+'R','${UV}'.split()))"`
UV_RE=`python -c "print ' '.join(map(lambda x: x+'E','${UV_R}'.split()[1:-1]))"`
UV_RE_R=`python -c "print ' '.join(map(lambda x: x[:-2],'${UV_RE}'.split()))"`

xrfi_simple.py -c 0_130,755_777,1540,1704,1827,1868,1885_2047 --df=6 $UV
# XXX is the "pol" argument below necessary?
ddr_filter_coarse.py -p xx --clean=1e-3 --maxbl=300 --output=ddr --invert $UV_R
xrfi_simple.py -n 4 $UV_RE
for FILE in $UV_RE_R; do
    echo Copying ${FILE}ER/flags to ${FILE}
    cp ${FILE}/flags ${FILE}/flags_bk
    cp ${FILE}ER/flags ${FILE}/flags
    #rm -rf ${FILE}[DEF]
done
rm -rf $UV_RE
ddr_filter_coarse.py -p xx --clean=1e-4 --maxbl=300 $UV_R
#rm -rf *.uvGCRER
#    #if ! ls ${FILE}drxam &> /dev/null; then
#    #fi
