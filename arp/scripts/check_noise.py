#! /usr/bin/env python
import numpy as np, pylab as plt
import capo, aipy
import sys, optparse

o = optparse.OptionParser()
o.set_usage('check_noise.py [options] *.uv')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o, cal=True, pol=True)
o.add_option('--skip_ants', dest='skip_ants', default='',
            help='Antennas to exclude from diagnostic.')
o.add_option('--verbose', action='store_true', help='Toggle verbosity')
opts,args = o.parse_args(sys.argv[1:])

aa = aipy.cal.get_aa(opts.cal, np.array([.15]))
capo.hex.aa_to_info_hera(aa)
info = capo.hex.aa_to_info_hera(aa)
reds = info.get_reds()
def get_red(sep):
     for r in reds:
          if sep in r: return r
all_bls = reduce(lambda x,y: x+y, reds)

POL = opts.pol
SKIP_ANTS = map(int, opts.skip_ants.split(','))
files = args
antstr = 'all'
uv = aipy.miriad.UV(files[0])
sdf, sfreq, nchan = uv['sdf'], uv['sfreq'], uv['nchan']
inttime = uv['inttime']
info, data, flgs = capo.miriad.read_files(files, antstr, POL)
for bl in data:
    data[bl] = data[bl][POL]
    flgs[bl] = flgs[bl][POL]
errs = capo.metrics.check_ants(reds, data, skip_ants=SKIP_ANTS, flag_thresh=.4)
bad_ants = [i for i in errs if errs[i] > len(errs)/2] + SKIP_ANTS
all_bls = [(i,j) for i,j in all_bls 
           if not i in bad_ants and not j in bad_ants]
print 'Bad antennas:', bad_ants

autos = [bl for bl in data.keys() 
         if bl[0] == bl[1] and not bl[0] in bad_ants]
print autos

xrfi = 0
for bl in autos:
    print 'Flagging', bl
    xrfi += capo.xrfi.xrfi(data[bl].real)
capo.plot.waterfall(xrfi, mode='lin', cmap='hot')
plt.xlabel('Channel')
plt.ylabel('Integration')
plt.title('RFI Flagging')
plt.show()
wgts = np.where(xrfi > 0, 0, 1.)

for i,j in all_bls:
    print 'Processing', i,j
    try: d = data[(i,j)]
    except(KeyError): d = data[(j,i)].conj()
    ai, aj = data[(i,i)].real, data[(j,j)].real
    ww = wgts[1:,1:] * wgts[1:,:-1] * wgts[:-1,1:] * wgts[:-1,:-1]
    dd = ((d[:-1,:-1] - d[:-1,1:]) - (d[1:,:-1] - d[1:,1:])) * ww / np.sqrt(4)
    dai = ((ai[:-1,:-1] + ai[:-1,1:]) + (ai[1:,:-1] + ai[1:,1:])) * ww / 4
    daj = ((aj[:-1,:-1] + aj[:-1,1:]) + (aj[1:,:-1] + aj[1:,1:])) * ww / 4
    Cij = np.sum(np.abs(dd)**2,axis=0)/np.sum(dai*daj,axis=0) * (.1e9/1024 * 10.737)
    plt.plot(Cij)
plt.ylim(0,5)
plt.grid()
plt.xlabel('Frequency [GHz]')
plt.ylabel('$|V_{N,ij}|^2/V_{ii}\cdot V_{jj}$')
plt.title('Normalized Noise')
plt.show()
