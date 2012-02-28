beamcal: input source track files, output beama (crosspoints beam) beamb (srctrack beam), srcflux v0 w/ lsq
mk_smooth_beam: input beam fits files, output smooth beam, smooth beam npz ("uv")
beam_renorm: input src track files, beam npz, output srctrack beam, srcflux v1

source tracks:
    <src_name>__<starttime>_<endtime>.<id>.npz
        id: s1cN = 1jy noise, f = real tracks (flagged)
    <src_name>.s__<starttime>_<endtime>.<id>.npz = uv-sim tracks

../beamcal_v12.py -C pgb015_v008 -p yy -o test.fits {cyg,cas,crab,vir}*.f.npz
> output test.fits.{npz,pkl}, test[ab].fits.fits
../mk_smooth_beam.py testb.fits.fits
> output output.npz

Simulations done with:
cyg is 10622.92
cas is 8914.3854135739166
crab is 1817.0860545717171
vir is 1400.7651852892755

Embedded xyz coordinates directly in npz files (originals in data/no_xyz)
