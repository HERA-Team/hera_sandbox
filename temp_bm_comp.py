#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import glob, pfits, ephem, sys

# Collect temperature data
flist = glob.glob('/data1/paper/gb/weather/2009_05_[12]*.fits')
flist.sort()
jds, temps = [], []
for f in flist:
    print f
    hdu = pfits.FITS(f).get_hdus()[1]
    dat = hdu.get_data()
    jds.append(dat['DMJD'] + 2400000.5)
    temps.append(dat['TEMPERATURE'])
jds = n.concatenate(jds)
temps = n.concatenate(temps)
JDBIN = ephem.second * 600
nbins = int((jds[-1]-jds[0])/JDBIN)
jwgts,bins = n.histogram(jds, bins=nbins)
jtemps,bins = n.histogram(jds, weights=temps, bins=nbins)
temp_prof = jtemps/jwgts

flist = sys.argv[1:]
flist.sort()

CH = 300
times, lsts, dat = [], [], []
for f in flist:
    print f
    uv = a.miriad.UV(f)
    #uv.select('antennae', 1, 3)
    uv.select('antennae', 1, 10)
    #uv.select('antennae', 1, 13)
    #uv.select('antennae', 0, 10)
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        times.append(t)
        lsts.append(uv['lst'])
        dat.append(d)
times = n.array(times)
lsts = n.array(lsts)
dat = n.abs(n.ma.array(dat))
dat = n.ma.masked_where(dat > 1, dat)

LSTBIN = .005
nbins = int(2*n.pi/LSTBIN)+1
def get_avg_prof(lsts, dat):
    avg_prof = []
    lwgts,lbins = n.histogram(lsts, range=(0,2*n.pi), bins=nbins)
    lwgts = lwgts.clip(1,n.inf)
    for ch in range(dat.shape[1]):
        d = dat[:,ch]
        nmsk = n.logical_not(d.mask)
        vlsts = lsts.compress(nmsk)
        vdat = d.compressed()
        lwgts,lbins = n.histogram(vlsts, range=(0,2*n.pi), bins=nbins)
        lwgts = lwgts.clip(1,n.inf)
        ldat,lbins = n.histogram(vlsts, weights=vdat, 
            range=(0,2*n.pi), bins=nbins)
        avg_prof.append(ldat/lwgts)
    return n.array(avg_prof).transpose()

# Get rid of some RFI
lstinds = n.floor(lsts / LSTBIN).astype(n.int)
sig = n.inf
for i in range(10):
    avg_prof = get_avg_prof(lsts, dat)
    avg_dat = avg_prof[lstinds,:]
    dif_dat = dat - avg_dat
    #if i == 0:
    #    p.subplot(131)
    #    p.imshow(n.log10(n.abs(dat)), vmax=0, vmin=-3, aspect='auto')
    #    p.colorbar(shrink=.5)
    #    p.subplot(132)
    #    p.imshow(n.log10(n.abs(dif_dat)), vmax=0, vmin=-3, aspect='auto')
    #    p.colorbar(shrink=.5)
    prv_sig = sig
    sig = n.ma.std(dif_dat)
    print sig, prv_sig
    if sig > prv_sig / 2: break
    dat = n.ma.masked_where(n.abs(dif_dat)/sig > 2, dat)

#p.subplot(133)
#p.imshow(n.log10(n.abs(dif_dat)), vmax=0, vmin=-3, aspect='auto')
#p.colorbar(shrink=.5)
#p.show()
#dat = dat.compress(valid)
#lsts = lsts.compress(valid)
#times = times.compress(valid)

lstinds = n.floor(lsts / LSTBIN).astype(n.int)
#lwgts,lbins = n.histogram(lsts, range=(0,2*n.pi), bins=nbins)
inds = n.floor((times - jds[0]) / JDBIN).astype(n.int)
temps = temp_prof[inds]
temps -= n.average(temps)

def gain(temps, H=-0.02): return 10**(H*temps/10)

def fitfunc(args):
    H = args[0]
    gains = gain(temps, H)
    chisq = 0
    for ch in range(dat.shape[1]):
        d = dat[:,ch] / gains
        nmsk = n.logical_not(d.mask)
        #nmsk = n.ones(d.shape)
        vlsts = lsts.compress(nmsk)
        if len(vlsts) == 0: continue
        vdat = d.compressed()
        lwgts,lbins = n.histogram(vlsts, range=(0,2*n.pi), bins=nbins)
        lwgts = n.where(lwgts == 0, 1, lwgts)
        ldat,lbins = n.histogram(vlsts, weights=vdat, 
            range=(0,2*n.pi), bins=nbins)
        avg_prof = ldat/lwgts
        avg_dat = avg_prof[lstinds]
        chisq += n.ma.average((d - avg_dat)**2)
        #print ch, chisq
        #if ch == 42:
        #    print lsts.shape, d.shape, vlsts.shape, vdat.shape
        #    print gains, d.shape, lwgts.shape
    chisq /= dat.shape[1]
    print H, chisq
    return chisq
    
H = a.optimize.fmin(fitfunc, [-0.01])
gains = gain(temps, H)
ldat,lbins = n.histogram(lsts, weights=dat/gains, 
    range=(0,2*n.pi), bins=nbins)
avg_prof = ldat/lwgts
p.subplot(211)
p.semilogy(times, gains, '^')
p.semilogy(times, dat, '.')
p.semilogy(times, dat/gains, '.')

p.subplot(212)
p.plot(lsts, dat, '.')
p.plot(lsts, dat/gains, '.')
p.plot((lbins[:-1] + lbins[1:])/2, avg_prof, '-')
p.show()
