"""
complist_gleam04.py

script for making gleam 045514.3 -300646
model component list in CASA
"""
import os
import shutil
import numpy as np
import argparse

args = argparse.ArgumentParser(description="Run with casa as: casa -c complist_gleam04.py <args>")
args.add_argument("-c", type=str, help="name of this script")
args.add_argument("--image", default=False, action='store_true', help='make FITS image of model')
args.add_argument("--freqs", default=None, type=str, help="comma-separated values for input into np.linspace({},{},{})")
args.add_argument("--cell", default='45arcsec', type=str, help="image pixel size in arcsec")
args.add_argument("--imsize", default=256, type=int, help="number of pixels in image")

if __name__ == "__main__":
    a = args.parse_args()

    # set variables
    direction = "J2000 04h55m14.3s -30d06m46s"
    ref_freq = "151MHz"

    # gleam 045514.3 -300646
    cl.addcomponent(label="GLEAM0455-3006", flux=17.1, fluxunit="Jy", dir=direction, freq=ref_freq,
                    shape="point", spectrumtype='spectral index', index=-0.78)

    # gleam 045826.5 -300717
    cl.addcomponent(label="GLEAM0458-3007", flux=12.7, fluxunit="Jy", dir="J2000 04h58m26.5 -30d07m17s", freq=ref_freq,
                    shape="point", spectrumtype='spectral index', index=-0.80)

    # gleam 044437.3 -280949
    cl.addcomponent(label="GLEAM0444-2809", flux=37.3, fluxunit="Jy", dir="J2000 04h44m37.3s -28d09m49s", freq=ref_freq,
                    shape="point", spectrumtype='spectral index', index=-0.76)


    # save
    if os.path.exists("gleam04.cl"):
        shutil.rmtree("gleam04.cl")
    cl.rename("gleam04.cl")

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
        ia.fromshape("gleam04.cl.im", [a.imsize, a.imsize, 1, Nfreqs], overwrite=True)
        cs = ia.coordsys()
        cs.setunits(['rad','rad','','Hz'])

        # set pixel properties
        cell_rad = qa.convert(qa.quantity(a.cell),"rad")['value']
        cs.setincrement([-cell_rad, cell_rad], type='direction')
        cs.setreferencevalue([qa.convert("04h55m14.3s",'rad')['value'], qa.convert("-30d06m46s",'rad')['value']], type="direction")

        # set freq properties
        qa_freqs = qa.quantity(freqs, 'MHz')
        cs.setspectral(frequencies=qa_freqs)
 
        # set flux properties, make image, export to fits
        ia.setcoordsys(cs.torecord())
        ia.setbrightnessunit("Jy/pixel")
        ia.modify(cl.torecord(), subtract=False)
        print("...saving gleam04.cl.fits")
        exportfits(imagename="gleam04.cl.im", fitsimage="gleam04.cl.fits", overwrite=True, stokeslast=False)

    cl.close()

