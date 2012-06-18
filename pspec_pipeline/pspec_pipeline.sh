#! /usr/bin/env bash

CAL=psa746_v010
RFI_CHAN=0_130,755_777,1540,1704,1827,1868,1901_2047
HORIZON=1.
for FILE in $*; do
    # from uvcbt
    xrfi_simple.py -a 0 --dt=4 --df=6 -c $RFI_CHAN --combine -t 20 $FILE
    pspec_prep.py -C $CAL --horizon=$HORIZON ${FILE}R
    xrfi_simple.py -n 5 ${FILE}RB
    pspec_to_npz_quick_n_dirty.py -p xx -C $CAL ${FILE}RBR
done
