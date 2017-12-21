"""
complist_BS.py

script for making Bernardi Sources
model component list in CASA
"""
import os
import shutil
import numpy as np
import argparse
import sys

args = argparse.ArgumentParser(description="Run with casa as: casa -c complist_BS.py <args>")
args.add_argument("-c", type=str, help="name of this script")
args.add_argument("--image", default=False, action='store_true', help='make FITS image of model')
args.add_argument("--freqs", default=None, type=str, help="comma-separated values for input into np.linspace({},{},{})")
args.add_argument("--cell", default='45arcsec', type=str, help="image pixel size in arcsec")
args.add_argument("--imsize", default=256, type=int, help="number of pixels in image")

if __name__ == "__main__":
    a = args.parse_args()

    # set variables
    direction = "J2000 21h01m38.3s -28d00m19s"
    ref_freq = "151MHz"

    # GLEAM 2101 -2800
    cl.addcomponent(label="GLEAM2101-2800", flux=16.9, fluxunit="Jy", dir=direction, freq=ref_freq,
                    shape="point", spectrumtype='spectral index', index=-0.76)

    # GLEAM 2107-2525
    cl.addcomponent(label="GLEAM2107-2525", flux=31.0, fluxunit="Jy", dir="J2000 21h07m22s -25d25m56s", freq=ref_freq,
                    shape="point", spectrumtype='spectral index', index=-0.72)

    # GLEAM 2107 -2529
    cl.addcomponent(label="GLEAM2107-2529", flux=14.1, fluxunit="Jy", dir="J2000 21h07m24s -25d29m53s", freq=ref_freq,
                    shape="point", spectrumtype='spectral index', index=)

    # GLEAM 2101 -2803
    cl.addcomponent(label="GLEAM2101-2803", flux=7.7, fluxunit="Jy", dir="J2000 21h01m41s -28d03m27s", freq=ref_freq,
                    shape="point", spectrumtype='spectral index', index=-0.72)

    # GLEAM 2100 -2829
    cl.addcomponent(label="GLEAM2100-2829", flux=4.1, fluxunit="Jy", dir="J2000 21h00m12s -28d29m04s", freq=ref_freq,
                    shape="point", spectrumtype='spectral index', index=-0.73)


    # save
    if os.path.exists("BS.cl"):
        shutil.rmtree("BS.cl")
    cl.rename("BS.cl")

    # make image
    if a.image:
        # get frequencies
        if a.freqs is None:
            Nfreqs = 1
            freqs = np.array([150.0])
        else:
            freqs = np.linspace(*np.array(a.freqs.split(',')).astype(np.float))
            Nfreqs = len(freqs)

        # setup image
        ia.fromshape("BS.cl.im", [a.imsize, a.imsize, 1, Nfreqs], overwrite=True)
        cs = ia.coordsys()
        cs.setunits(['rad','rad','','Hz'])

        # set pixel properties
        cell_rad = qa.convert(qa.quantity(a.cell),"rad")['value']
        cs.setincrement([-cell_rad, cell_rad], type='direction')
        cs.setreferencevalue([qa.convert("21h01m38.3s",'rad')['value'], qa.convert("-28d00m19s",'rad')['value']], type="direction")

        # set freq properties
        qa_freqs = qa.quantity(freqs, 'MHz')
        cs.setspectral(frequencies=qa_freqs)
 
        # set flux properties, make image, export to fits
        ia.setcoordsys(cs.torecord())
        ia.setbrightnessunit("Jy/pixel")
        ia.modify(cl.torecord(), subtract=False)
        exportfits(imagename="BS.cl.im", fitsimage="BS.cl.fits", overwrite=True, stokeslast=False)

    cl.close()





