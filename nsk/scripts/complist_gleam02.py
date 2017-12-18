"""
complist_gleam02.py

script for making gleam 020012 -305327
model component list in CASA
"""
import os
import shutil
import numpy as np
import argparse

args = argparse.ArgumentParser(description="Run with casa as: casa -c complist_gleam02.py <args>")
args.add_argument("-c", type=str, help="name of this script")
args.add_argument("--image", default=False, action='store_true', help='make FITS image of model')

if __name__ == "__main__":
    a = args.parse_args()

    # set variables
    direction = "J2000 02h00m12.7s -30d53m27s"
    ref_freq = "151MHz"

    # gleam 020012 -305327
    cl.addcomponent(label="GLEAM0200-3053", flux=17.5, fluxunit="Jy", dir=direction, freq=ref_freq,
                    shape="point", spectrumtype='spectral index', index=-0.86)

    # gleam 0150 -2931
    cl.addcomponent(label="GLEAM0200-3053", flux=16.6, fluxunit="Jy", dir="J2000 01h50m36s -29d31m59s", freq=ref_freq,
                    shape="point", spectrumtype='spectral index', index=-0.80)

    # gleam 015200 -294056
    cl.addcomponent(label="GLEAM0200-3053", flux=4.7, fluxunit="Jy", dir="J2000 01h52m00s -29d40m56s", freq=ref_freq,
                    shape="point", spectrumtype='spectral index', index=-0.73)

    # save
    if os.path.exists("gleam02.cl"):
        shutil.rmtree("gleam02.cl")
    cl.rename("gleam02.cl")

    # make image
    if a.image:
        ia.fromshape("gleam02.cl.im", [512, 512, 1, 1], overwrite=True)
        cs = ia.coordsys()
        cs.setunits(['rad','rad','','Hz'])
        cell_rad = qa.convert(qa.quantity("45arcsec"),"rad")['value']
        cs.setincrement([-cell_rad, cell_rad], 'direction')
        cs.setreferencevalue([qa.convert("02h00m12.7s",'rad')['value'], qa.convert("-30d53m27s",'rad')['value']], type="direction")
        cs.setreferencevalue('151MHz','spectral')
        cs.setincrement('1MHz','spectral')
        ia.setcoordsys(cs.torecord())
        ia.setbrightnessunit("Jy/pixel")
        ia.modify(cl.torecord(), subtract=False)
        exportfits(imagename="gleam02.cl.im", fitsimage="gleam02.cl.fits", overwrite=True)

    cl.close()





