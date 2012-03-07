#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys

def get_dict_of_uv_data(filenames, antstr, polstr, decimate=1, decphs=0, verbose=False):
    times, dat, flg = [], {}, {}
    if type(filenames) == 'str': filenames = [filenames]
    for filename in filenames:
        if verbose: print '   Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, antstr, polstr)
        if decimate > 1: uv.select('decimate', decimate, decphs)
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            if len(times) == 0 or t != times[-1]:
                times.append(t)
            bl = a.miriad.ij2bl(i,j)
            dat[bl] = dat.get(bl,[]) + [d]
            flg[bl] = flg.get(bl,[]) + [f]
    for bl in dat:
        dat[bl] = n.array(dat[bl])
        flg[bl] = n.array(flg[bl])
    return n.array(times), dat, flg

rewire = {
     1: 16,  2: 17,  3: 18,  4: 19,
     5: 20,  6: 21,  7: 22,  8: 23,
     9:  0, 10:  1, 11:  2, 12:  3,
    13:  4, 14:  5, 15:  6, 16:  7,
    17:  8, 18:  9, 19: 10, 20: 11,
    21: 12, 22: 13, 23: 14, 24: 15,
}

INT = 10

times,dat,flg = get_dict_of_uv_data([sys.argv[-1]], 'all,-24,-25,-26,-27,-28,-29,-30,-31', 'yy')
NCHAN = dat.values()[0].shape[1]
print NCHAN

dout = []

for i in range(1,24):
    for j in range(i,24):
        _i = rewire[i]
        _j = rewire[j]
        bl = a.miriad.ij2bl(_i,_j)
        try: t,d,f = times[INT], dat[bl][INT], flg[bl][INT]
        except(KeyError): d = n.zeros((NCHAN,), dtype=n.complex64)
        dout.append(d)
        print len(dout)
dout = n.transpose(n.array(dout))
print dout.shape

for i in range(dout.shape[0]):
    for j in range(dout.shape[1]):
        print dout[i,j].real, dout[i,j].imag,
    print

        
