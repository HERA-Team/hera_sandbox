#! /usr/bin/env python
import aipy as a 
import numpy as n

aa = a.cal.get_aa('psa6240_v003', 0.0004926108374384237, .1, 203)
freqs = aa.get_afreqs()

h = a.healpix.HealpixMap(nside=64)
#get xyz cords of all patches.
xyz = h.px2crd(n.arange(h.map.size), ncrd=3)
top_x, top_y, top_z = n.dot(aa._eq2zen, xyz)
_bmx = aa[0].bm_response((top_x,top_y,top_z), pol='x')
_bmy = aa[0].bm_response((top_x,top_y,top_z), pol='y')
above_zero = n.where(top_z > 0)
bmxx = _bmx[:,above_zero].reshape(len(freqs), len(above_zero[0]))
bmyy = _bmy[:,above_zero].reshape(len(freqs), len(above_zero[0]))
top_x = top_x[above_zero]
top_y = top_y[above_zero]
top_z = top_z[above_zero]
az,alt = a.coord.top2azalt((top_x,top_y,top_z))





n.savez('PAPER_beam', bmxx=bmxx, bmyy=bmyy, freqs=freqs, az=az,alt=alt)
