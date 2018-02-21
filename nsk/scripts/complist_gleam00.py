"""
complist_gleam00.py

script for making gleam 002430.2 -292847
model component list in CASA
"""
import os
import shutil
import numpy as np
import argparse
import sys

args = argparse.ArgumentParser(description="Run with casa as: casa -c complist_gleam02.py <args>")
args.add_argument("-c", type=str, help="name of this script")
args.add_argument("--image", default=False, action='store_true', help='make FITS image of model')
args.add_argument("--freqs", default=None, type=str, help="comma-separated values for input into np.linspace({},{},{})")
args.add_argument("--cell", default='45arcsec', type=str, help="image pixel size in arcsec")
args.add_argument("--imsize", default=256, type=int, help="number of pixels in image")

if __name__ == "__main__":
    a = args.parse_args()

    # set variables
    direction = "J2000 00h24m30.2s -29d28m47s"
    ref_freq = "151MHz"

    # gleam 0024 -2928
    cl.addcomponent(label="GLEAM0024-2928", flux=16.1, fluxunit="Jy", dir=direction, freq=ref_freq,
                    shape="point", spectrumtype='spectral index', index=-0.86)

    # gleam 0025 -2602
    cl.addcomponent(label="GLEAM0025-2602", flux=17.6, fluxunit="Jy", dir="00h25m49s -26d02m11s",
                    freq=ref_freq, shape="point", spectrumtype='spectral index', index=-0.80)

    # gleam 0025 -3303
    cl.addcomponent(label="GLEAM0025-3303", flux=8.8, fluxunit="Jy", dir="00h25m30.5s -33d03m36s",
                    freq=ref_freq, shape="point", spectrumtype='spectral index', index=-0.86)

    # save
    if os.path.exists("gleam00.cl"):
        shutil.rmtree("gleam00.cl")
    cl.rename("gleam00.cl")

    # make image
    if a.image:
        # get frequencies
        if a.freqs is None:
            Nfreqs = 1
            freqs = np.array([151.0])
        else:
            freqs = np.linspace(*np.array(a.freqs.split(',')).astype(np.float))
            Nfreqs = len(freqs)

        # setup image
        ia.fromshape("gleam00.cl.im", [a.imsize, a.imsize, 1, Nfreqs], overwrite=True)
        cs = ia.coordsys()
        cs.setunits(['rad','rad','','Hz'])

        # set pixel properties
        cell_rad = qa.convert(qa.quantity(a.cell),"rad")['value']
        cs.setincrement([-cell_rad, cell_rad], type='direction')
        cs.setreferencevalue([qa.convert("00h24m30.2s",'rad')['value'], qa.convert("-29d28m47s",'rad')['value']], type="direction")

        # set freq properties
        qa_freqs = qa.quantity(freqs, 'MHz')
        cs.setspectral(frequencies=qa_freqs)
 
        # set flux properties, make image, export to fits
        ia.setcoordsys(cs.torecord())
        ia.setbrightnessunit("Jy/pixel")
        ia.modify(cl.torecord(), subtract=False)
        print("...saving gleam00.cl.fits")
        exportfits(imagename="gleam00.cl.im", fitsimage="gleam00.cl.fits", overwrite=True, stokeslast=False)

    cl.close()





