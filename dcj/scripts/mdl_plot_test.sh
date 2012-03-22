#!/bin/bash
#git stash apply
BUILD=build/scripts-2.6
python setup.py clean -a
python setup.py install --prefix=~
mdlvis.py --sim -C pgb966_v003 -s her --startjd=2453887 --endjd=2453887.02 --pol=yy --sfreq=0.156 --sdf=0.000234 --nchan=1 --inttime=60 & plot_uv.py -a 0,1 -p yy new.uv -o test.png
rm -r new.uv
#rv=$?
#
#exit $rv