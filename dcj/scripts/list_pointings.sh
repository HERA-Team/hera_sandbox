#! /bin/bash
#list files with their pointings, sorted by RA
#usage list_pointings.sh *c60*bim.fits 
find_pointing.py -v $* | sort -tJ -n -k2 > pointing_list.txt