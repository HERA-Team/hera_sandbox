#! /bin/bash 
CAT=vlss
WIDTH=1
extract_subimg.py -s all --cat=${CAT} $* --prefix=${CAT} --width=${WIDTH}
montage +frame +shadow +label -geometry 200x200+0+0 ${CAT}*_thumb.fits -background grey25 ${CAT}_montage.png