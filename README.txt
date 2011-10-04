beamcal: input source track files, output beama (crosspoints beam) beamb (srctrack beam), srcflux v0 w/ lsq
mk_smooth_beam: input beam fits files, output smooth beam, smooth beam npz ("uv")
beam_renorm: input src track files, beam npz, output srctrack beam, srcflux v1

source tracks:
    <src_name>__<starttime>_<endtime>.<id>.npz
        id: s1cN = 1jy noise, f = real tracks (flagged)
    <src_name>.s__<starttime>_<endtime>.<id>.npz = uv-sim tracks
