#! /bin/bash

# Pre-lstbin pipeline
xrfi_simple.py -c 0_139,1900_2047 --combine --thresh=10 --df=6 --dt=4 $*

# Post-lstbin pipelin
for file in $* ; do
    ~/pkgs/capo/pspec_pipline/pspec_prep.py -C psa898_v002 --nogain --nophs --window=none --horizon=2.0 ${file}
    xtalk3.py ${file}B
    xrfi_simple.py -c 0_170,341,647_649,742,744_746,1076,1077,1537,1538,1662_1664,1685,1701_1707,1737,1767,1795,1799,1800,1802,1807,1808,1825,1847,1848,1855_2047 --dt=4 ${file}Bx
    ~/pkgs/capo/pspec_pipline/pspec_prep.py -C psa898_v002 --nogain --horizon=2.0 ${file}BxR
    xrfi_simple.py -n 4 ${file}BxRB
    ~/pkgs/capo/pspec_pipline/pspec_to_npz_quick_n_dirty.py -C psa898_v002 -p xx ${file}BxRBR
done


    
