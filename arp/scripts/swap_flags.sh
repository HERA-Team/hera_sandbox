#! /bin/bash

for FILE in $* ; do
    xrfi_simple.py -c 0_130,755_777,1540,1704,1827,1868,1901_2047 --dt=4 --df=6 ${FILE}
    ~/pkgs/capo/pspec_pipline/pspec_prep.py -C psa898_v002 --nogain --nophs --horizon=2.0 ${FILE}R
    xtalk3.py ${FILE}RB
    xrfi_simple.py -n 4 --combine --thresh=1 ${FILE}RBx
    mv ${FILE}R/flags ${FILE}R/flags_bk
    cp ${FILE}RBxR/flags ${FILE}R/flags
    #mv $FILE/flags_bk $FILE/flags
done
