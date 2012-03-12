#! /bin/bash
#Usage make_cube.sh *.fits
#Input a set of channel maps and it merges them and magically fixes the headers
fix_arp_fits.py $*
FIXED=`fix_arp_fits.py -l $*`
stack_fits.py $FIXED
STACK=`find_pointing.py -J newstack.fits`_cube.fits
echo creating $STACK
mv newstack.fits $STACK
cube_crop.py $STACK

