"""
sky_reduce.py
-------------

run with casa as:
casa -c sky_reduce.py <args>
"""
import sys
import os
import numpy as np
import argparse
import subprocess
import shutil
import glob

# get args
a = argparse.ArgumentParser()
a.add_argument('--script', '-c', type=str, help='name of this script')
a.add_argument('--msin', default=None, type=str, help='path to CASA measurement set. if fed a .uvfits, will convert to ms', required=True)
a.add_argument('--source', default=None, type=str, help='source name', required=True)
a.add_argument('--refant', default=None, type=str, help='reference antenna', required=True)
a.add_argument('--ex_ants', default=None, type=str, help='bad antennas to flag')
a.add_argument('--antmets', default=None, type=str, help='path to antenna metrics file to use in antenna flagging')
a.add_argument('--unflag', default=False, action='store_true', help='start by unflagging data')
a.add_argument('--rflag', default=False, action='store_true', help='run flagdata(mode=rflag)')
a.add_argument('--out_dir', default=None, type=str, help='output directory')
a.add_argument('--nocal', default=False, action='store_true', help='skip calibration and just make an image')
args = a.parse_args()

# get vars
source = args.source
refant = args.refant
ex_ants = args.ex_ants
rflag = args.rflag
antmets = args.antmets
nocal = args.nocal

# get paths
msin = args.msin
base_ms = os.path.basename(msin)
if args.out_dir is None:
    out_dir = os.path.dirname(msin)
else:
    out_dir = args.out_dir

# check for uvfits
if base_ms.split('.')[-1] == 'uvfits':
    print("...converting uvfits to ms")
    uvfits = msin
    msin = os.path.join(out_dir, '.'.join(base_ms.split('.')[:-1] + ['ms']))
    base_ms = os.path.basename(msin)
    if os.path.exists(msin):
        shutil.rmtree(msin)
    importuvfits(uvfits, msin)

# fix vis
print("...fix vis")
ra, dec = np.loadtxt('{}.loc'.format(source), dtype=str)
ra, dec = ra.split(':'), dec.split(':')
fixvis(msin, msin, phasecenter='J2000 {}h{}m{}s {}d{}m{}s'.format(*(ra+dec)))

# insert source model
# set up calibrator source
try:
    cl.addcomponent(flux=1, fluxunit='Jy', shape='point', dir='J2000 {}h{}m{}s {}d{}m{}s'.format(*(ra+dec)))
    cl.rename(os.path.join(out_dir, '{}.cl'.format(source)))
    cl.close()
except:
    pass

ft(msin, complist=os.path.join(out_dir, "{}.cl".format(source)), usescratch=True)

# unflag
if args.unflag is True:
    print("...unflagging")
    flagdata(msin, mode='unflag')

# flag bad ants
if ex_ants is not None:
    print("...flagging bad ants: {}".format(ex_ants))
    flagdata(msin, mode='manual', antenna=ex_ants)

# rflag
if rflag is True:
    print("...rfi flagging")
    flagdata(msin, autocorr=True)
    flagdata(msin, mode='rflag')

# setup calibration tables
kc = os.path.join(out_dir, base_ms+'.{}.cal'.format('K'))
gc = os.path.join(out_dir, base_ms+'.{}.cal'.format('G'))

if nocal is True:
    Bsplit = msin
else:
    # perform initial K calibration (per-antenna delay)
    print("...performing K gaincal")
    if os.path.exists(kc):
        shutil.rmtree(kc)
    if os.path.exists("{}.png".format(kc)):
        os.remove("{}.png".format(kc))
    gaincal(msin, caltable=kc, gaintype="K", solint='inf', refant=refant, minsnr=1, spw="0:150~900")
    plotcal(kc, xaxis='antenna', yaxis='delay', figfile='{}.png'.format(kc), showgui=False)

    # write delays to text file
    tb.open(kc)
    delays = tb.getcol('FPARAM')[0, 0][~tb.getcol('FLAG')[0, 0]]
    delay_ants = tb.getcol('ANTENNA1')[~tb.getcol('FLAG')[0, 0]]
    tb.close()
    np.savetxt("{}.csv".format(kc), np.vstack([delay_ants, delays]).T, fmt="%12.8f", header="ant_num, delay [ns]", delimiter=", ")

    # perform initial G calibration (per-spw and per-pol gain)
    print("...performing G gaincal")
    if os.path.exists(gc):
        shutil.rmtree(gc)
    gaincal(msin, caltable=gc, gaintype='G', solint='inf', refant=refant, minsnr=2, calmode='p', spw="0:150~900", gaintable=kc)

    # apply calibrations
    print("...applying KG gaincal")
    applycal(msin, gaintable=[kc, gc])

    # split data
    KGsplit = os.path.join(out_dir, "{}.KGsplit".format(base_ms))
    files = glob.glob("{}*".format(KGsplit))
    for f in files:
        if os.path.exists(f):
            try:
                shutil.rmtree(f)
            except OSError:
                os.remove(f)
    split(msin, KGsplit, datacolumn="corrected")

    # calibrated bandpass
    print("...performing B bandpass cal")
    bc = KGsplit+'.{}.cal'.format('B')
    bandpass(vis=KGsplit, spw="", minsnr=1, solnorm=False, bandtype='B', caltable=bc)
    plotcal(bc, xaxis='chan', yaxis='amp', figfile="{}.amp.png".format(bc), showgui=False)
    plotcal(bc, xaxis='chan', yaxis='phase', figfile="{}.phs.png".format(bc), showgui=False)

    applycal(KGsplit, gaintable=[bc])
    Bsplit = os.path.join(out_dir, "{}.Bsplit".format(base_ms))
    files = glob.glob("{}*".format(Bsplit))
    for f in files:
        if os.path.exists(f):
            try:
                shutil.rmtree(f)
            except OSError:
                os.remove(f)

    split(KGsplit, Bsplit, datacolumn="corrected")

# clean
print("...performing clean")
im_stem = os.path.join(out_dir, base_ms + '.' + source)
source_files = glob.glob(im_stem+'*')
if len(source_files) > 0:
    for sf in source_files:
        if os.path.exists(sf):
            try:
                shutil.rmtree(sf)
            except OSError:
                os.remove(sf)

clean(vis=Bsplit, imagename=im_stem, spw="0:200~850", niter=500, weighting='briggs', robust=-0.5, imsize=[512, 512],
      cell=['300arcsec'], mode='mfs', nterms=1)#, mask='circle[[{}h{}m{}s, {}d{}m{}s ], 7deg]'.format(*(ra+dec)))
exportfits(imagename='{}.image'.format(im_stem), fitsimage='{}.fits'.format(im_stem))

