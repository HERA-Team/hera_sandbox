"""
complist_gleam05.py

script for making gleam 052257.9 -362727
model component list in CASA
"""
import os
import shutil
import numpy as np
import argparse

args = argparse.ArgumentParser(description="Run with casa as: casa -c complist_gleam05.py <args>")
args.add_argument("-c", type=str, help="name of this script")
args.add_argument("--image", default=False, action='store_true', help='make FITS image of model')
args.add_argument("--freqs", default=None, type=str, help="comma-separated values for input into np.linspace({},{},{})")
args.add_argument("--cell", default='45arcsec', type=str, help="image pixel size in arcsec")
args.add_argument("--imsize", default=256, type=int, help="number of pixels in image")

if __name__ == "__main__":
    a = args.parse_args()

    # set variables
    direction = "J2000 05h22m57.9s -36d27m27s"
    ref_freq = "151MHz"

    # gleam 052257.9 -362727
    cl.addcomponent(label="GLEAM0522-3627", flux=55.9, fluxunit="Jy", dir=direction, freq=ref_freq,
                    shape="point", spectrumtype='spectral index', index=-0.62)

    # gleam 053949.1 -341235
    cl.addcomponent(label="GLEAM0539 -3412", flux=16.6, fluxunit="Jy", dir="J2000 05h39m49.1s -34d12m35s", freq=ref_freq,
                    shape="point", spectrumtype='spectral index', index=-0.94)


    # save
    if os.path.exists("gleam05.cl"):
        shutil.rmtree("gleam05.cl")
    cl.rename("gleam05.cl")

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
        ia.fromshape("gleam05.cl.im", [a.imsize, a.imsize, 1, Nfreqs], overwrite=True)
        cs = ia.coordsys()
        cs.setunits(['rad','rad','','Hz'])

        # set pixel properties
        cell_rad = qa.convert(qa.quantity(a.cell),"rad")['value']
        cs.setincrement([-cell_rad, cell_rad], type='direction')
        cs.setreferencevalue([qa.convert("05h22m57.9s",'rad')['value'], qa.convert("-36d27m27s",'rad')['value']], type="direction")

        # set freq properties
        qa_freqs = qa.quantity(freqs, 'MHz')
        cs.setspectral(frequencies=qa_freqs)
 
        # set flux properties, make image, export to fits
        ia.setcoordsys(cs.torecord())
        ia.setbrightnessunit("Jy/pixel")
        ia.modify(cl.torecord(), subtract=False)
        print("...saving gleam05.cl.fits")
        exportfits(imagename="gleam05.cl.im", fitsimage="gleam05.cl.fits", overwrite=True, stokeslast=False)

    cl.close()

