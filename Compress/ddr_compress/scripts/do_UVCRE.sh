#! /bin/bash

rm -rf $2E
ddr_filter_coarse.py -a 1 -p xx,xy,yx,yy --clean=1e-3 --maxbl=300 --output=ddr --invert $1 $2 $3
