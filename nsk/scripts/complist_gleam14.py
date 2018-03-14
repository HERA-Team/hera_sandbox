"""
complist_gleam14.py

script for making gleam 14h25m28.9s -29d59m55s
model component list in CASA
"""
import os
import shutil
import numpy as np
import argparse
import sys

args = argparse.ArgumentParser(description="Run with casa as: casa -c complist_gleam14.py <args>")
args.add_argument("-c", type=str, help="name of this script")
args.add_argument("--image", default=False, action='store_true', help='make FITS image of model')
args.add_argument("--freqs", default=None, type=str, help="comma-separated values for input into np.linspace({},{},{})")
args.add_argument("--cell", default='45arcsec', type=str, help="image pixel size in arcsec")
args.add_argument("--imsize", default=256, type=int, help="number of pixels in image")

if __name__ == "__main__":
    a = args.parse_args()

    # set variables
    direction = "J2000 14h25m28.9s -29d59m55s"
    ref_freq = "151MHz"

    # gleam 1425 -2959
    cl.addcomponent(label="GLEAM1425-2959", flux=16.0, fluxunit="Jy", dir=direction, freq=ref_freq,
                    shape="point", spectrumtype='spectral index', index=-0.90)

    # gleam 1421 -3104
    cl.addcomponent(label="GLEAM1421-3104", flux=7.8, fluxunit="Jy", dir="14h21m55.4s -31d04m22s",
                    freq=ref_freq, shape="point", spectrumtype='spectral index', index=-0.85)

    # gleam 1422 -2727
    cl.addcomponent(label="GLEAM1422-2727", flux=19.7, fluxunit="Jy", dir="14h22m49.4s -27d27m55s",
                    freq=ref_freq, shape="point", spectrumtype='spectral index', index=-0.99)

    # gleam 1442 -2637
    cl.addcomponent(label="GLEAM1442-2637", flux=9.6, fluxunit="Jy", dir="14h42m02s -26d37m12s",
                    freq=ref_freq, shape="point", spectrumtype='spectral index', index=-0.84)

    # gleam 1421 -3104
    cl.addcomponent(label="GLEAM1421-3104", flux=7.8, fluxunit="Jy", dir="14h21m55s -31d04m22s",
                    freq=ref_freq, shape="point", spectrumtype='spectral index', index=-0.9)

    # save
    if os.path.exists("gleam14.cl"):
        shutil.rmtree("gleam14.cl")
    cl.rename("gleam14.cl")

    # make image
    if a.image:
        # get frequencies
        if a.freqs is None:
            Nfreqs = 1
            freqs = np.array([151.0])
        else:
            freqs = np.linspace(*np.array(a.freqs.split(',')).astype(np.float), endpoint=True)
            Nfreqs = len(freqs)

        # setup image
        ia.fromshape("gleam14.cl.im", [a.imsize, a.imsize, 1, Nfreqs], overwrite=True)
        cs = ia.coordsys()
        cs.setunits(['rad','rad','','Hz'])

        # set pixel properties
        cell_rad = qa.convert(qa.quantity(a.cell),"rad")['value']
        cs.setincrement([-cell_rad, cell_rad], type='direction')
        cs.setreferencevalue([qa.convert("14h25m28.9s",'rad')['value'], qa.convert("-29d59m55s",'rad')['value']], type="direction")

        # set freq properties
        qa_freqs = qa.quantity(freqs, 'MHz')
        cs.setspectral(frequencies=qa_freqs)
 
        # set flux properties, make image, export to fits
        ia.setcoordsys(cs.torecord())
        ia.setbrightnessunit("Jy/pixel")
        ia.modify(cl.torecord(), subtract=False)
        print("...saving gleam14.cl.fits")
        exportfits(imagename="gleam14.cl.im", fitsimage="gleam14.cl.fits", overwrite=True, stokeslast=False)

    cl.close()





