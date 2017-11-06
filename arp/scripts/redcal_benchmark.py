#! /usr/bin/env python

import hera_cal.redcal as redcal
import omnical
import numpy as np
import cProfile

np.random.seed(0)

def build_reds_hex(hexNum, sep=14.7):
    antpos, i = {}, 0
    for row in range(hexNum-1,-(hexNum),-1):
        for col in range(2*hexNum-abs(row)-1):
            xPos = ((-(2*hexNum-abs(row))+2)/2.0 + col)*sep;
            yPos = row*sep*3**.5/2;
            antpos[i] = np.array([xPos, yPos, 0])
            i += 1
    return redcal.get_reds(antpos), antpos

class RedundantInfo(omnical.calib.RedundantInfo):
    def order_data(self, dd):
        d = []
        for i, j in self.bl_order():
            bl = (i, j)
            pol = 'xx'
            try:
                d.append(dd[bl+(pol,)])
            except(KeyError):
                d.append(dd[bl[::-1]+(pol[::-1],)].conj())
        rv = np.array(d).transpose((1, 2, 0))
        return rv

# Overhead
RCOND = 1e-6
HEXNUM = 3
reds, antpos = build_reds_hex(HEXNUM)
nants = len(antpos)
antpos_omni = np.empty((nants,3), dtype=np.float)
for i in antpos: antpos_omni[i] = antpos[i]
print '# Antennas:', nants

dtype = np.complex64
#dtype = np.complex128
gains, true_vis, d = redcal.sim_red_data(reds, ['xx'], shape=(30,30), gain_scatter=.01)
d = {k:dk.astype(dtype) for k,dk in d.items()}
w = {k:1. for k in d.keys()}
sol0 = {k:np.ones(v.shape, dtype=dtype) for k,v in gains.items()}


rc = redcal.RedundantCalibrator(reds, antpos)
sol0.update(rc.compute_ubls(d,sol0))
info = RedundantInfo()
info.init_from_reds(reds, antpos_omni)

cProfile.run('meta1a, sol1a = rc.lincal(d, sol0, conv_crit=RCOND, sparse=False, verbose=False)', 'redcal_hex%d_redcal.prf' % HEXNUM)
cProfile.run('meta2b, gains1b, vis1b = omnical.calib.lincal(d, info, gains=sol0)', 'redcal_hex%d_omnical.prf' % HEXNUM)


print 'Done'
