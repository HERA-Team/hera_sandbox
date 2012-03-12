#! /bin/bash
NCHAN=256
XRFI_CH=0_59,780_1023

# ---------------------------------------

for FILE in $*; do
    echo Working on $FILE
    JD=`echo $FILE | python -c "import sys; print '.'.join(sys.stdin.read().split('.')[1:3])"`
    UVCB=zen.${JD}.uvcb
    if ! ls ${UVCB}rm &> /dev/null; then
        echo "File ${UVCB}rm doesn't exist"
        echo "Creating it..."
        xrfi.py -m val -c $XRFI_CH ${UVCB}
        #xtalk3.py ${UVCB}r
        combine_freqs.py -n $NCHAN ${UVCB}r
        #rm -r ${UVCB}r
    fi

done
