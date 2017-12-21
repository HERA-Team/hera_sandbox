"""
complist_fornaxA.py

script for making fornaxA
model component list in CASA
"""
import os
import shutil
import numpy as np
import argparse
import sys

args = argparse.ArgumentParser(description="Run with casa as: casa -c complist_fornaxA.py <args>")
args.add_argument("-c", type=str, help="name of this script")
args.add_argument("--image", default=False, action='store_true', help='make FITS image of model')
args.add_argument("--freqs", default=None, type=str, help="comma-separated values for input into np.linspace({},{},{})")
args.add_argument("--cell", default='45arcsec', type=str, help="image pixel size in arcsec")
args.add_argument("--imsize", default=256, type=int, help="number of pixels in image")

if __name__ == "__main__":
    a = args.parse_args()

    # set variables
    core_direction = "J2000 3h22m45s -37d12m00s"
    ref_freq = "154MHz"

    # West Lobe: McKinley et al. 2015
    cl.addcomponent(label="West Lobe", flux=480, fluxunit="Jy", dir="J2000 3h21m30s -37d10m00s", freq=ref_freq,
                    shape="gaussian", majoraxis="20arcmin", minoraxis="20arcmin", positionangle="0deg",
                    spectrumtype='spectral index', index=-0.77)

    # East Lobe: McKinley et al. 2015
    cl.addcomponent(label="East Lobe", flux=260, fluxunit="Jy", dir="J2000 3h24m00s -37d16m00s", freq=ref_freq,
                    shape="gaussian", majoraxis="15arcmin", minoraxis="15arcmin", positionangle="0deg",
                    spectrumtype='spectral index', index=-0.77)

    # Core: McKinley et al. 2015
    cl.addcomponent(label="Core", flux=12, fluxunit="Jy", dir=core_direction, freq=ref_freq,
                    shape="gaussian", majoraxis="5arcmin", minoraxis="5arcmin", positionangle="0deg",
                    spectrumtype='spectral index', index=-1.0)

    # save
    if os.path.exists("fornaxA.cl"):
        shutil.rmtree("fornaxA.cl")
    cl.rename("fornaxA.cl")

    # make image
    if a.image:
        # get frequencies
        if a.freqs is None:
            Nfreqs = 1
            freqs = np.array([154.0])
        else:
            freqs = np.linspace(*np.array(a.freqs.split(',')).astype(np.float))
            Nfreqs = len(freqs)

        # setup image
        ia.close()
        ia.fromshape("fornaxA.cl.im", [a.imsize, a.imsize, 1, Nfreqs], overwrite=True)
        cs = ia.coordsys()
        cs.setunits(['rad','rad','','Hz'])

        # set pixel properties
        cell_rad = qa.convert(qa.quantity(a.cell),"rad")['value']
        cs.setincrement([-cell_rad, cell_rad], type='direction')
        cs.setreferencevalue([qa.convert("03h22m45s",'rad')['value'], qa.convert("-37d12m00s",'rad')['value']], type="direction")

        # set freq properties
        qa_freqs = qa.quantity(freqs, 'MHz')
        cs.setspectral(frequencies=qa_freqs)
 
        # set flux properties, make image, export to fits
        ia.setcoordsys(cs.torecord())
        ia.setbrightnessunit("Jy/pixel")
        ia.modify(cl.torecord(), subtract=False)
        exportfits(imagename="fornaxA.cl.im", fitsimage="fornaxA.cl.fits", overwrite=True, stokeslast=False)

    cl.close()





