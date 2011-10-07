#! /bin/bash

echo Crossing Points
./beamcal.py $* -o test | grep flux
./mk_smooth_beam.py testb.fits > /dev/null
echo Sourcetrack renorm1
./beam_renorm.py -o testc -b testb_sm.npz $* | grep flux
./mk_smooth_beam.py testc.fits > /dev/null
echo Sourcetrack renorm2
./beam_renorm.py -o testd -b testc_sm.npz $* | grep flux
./mk_smooth_beam.py testd.fits > /dev/null
echo Sourcetrack renorm3
./beam_renorm.py -o teste -b testd_sm.npz $* | grep flux
rm teste.fits

