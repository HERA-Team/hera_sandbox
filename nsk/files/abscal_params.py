"""
abscal_params.py

abscal parameter file
"""
## set flags ##
T = True
F = False

overwrite           = T

run_abscal          = F
rfi_flag            = F
casa_rflag          = F

apply_abscal        = F
multiprocess        = F
Nproc               = 16
abs_xx_calfits      = "../data/H1C/2458042/zen.2458042.45361.xx.HH.uvR.ms.abs.calfits"
abs_yy_calfits      = "../data/H1C/2458042/zen.2458042.45361.yy.HH.uvR.ms.abs.calfits"

source_spectrum     = T
spec_rfi_flag       = F
spec_casa_rflag     = T
spec_dchan          = 5
gf_mult             = 1
rms_maxr            = 4
rms_minr            = 1

## set field ##
field = 'gleam04'
source = 'gleam0444'

## set calibration parameters ##
ex_ants = '0,2,11,14,50,98'
refant = '53'
duration = 4.0
imsize = 512
pxsize = 250
niter = 50

JD = 2458042
data_path = os.path.join("/lustre/aoc/projects/hera/nkern/data/H1C", str(JD))
complist = "complist_{}.py".format(field)
beamfile = "../Software/HERA-Beams/NicolasFagnoniBeams/NF_HERA_Beams.beamfits"
lon = 21.428303
lat = -30.721526

