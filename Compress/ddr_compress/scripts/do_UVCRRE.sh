#! /bin/bash

rm -rf $2E
ddr_filter_coarse.py -a all -p xx,xy,yx,yy --maxbl=301 --clean=1e-3 --nsections=20 $1 $2 $3
