#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import glob, pfits, ephem, sys

'''
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
Ts = jtemps/jwgts.clip(1,n.inf)
Ts -= n.average(Ts)

def gain(jd, H=-0.06): 
    ind = n.floor((jd - jds[0]) / JDBIN).astype(n.int)
    T = Ts[ind]
    return 10**(H*T/10)
'''

flist = sys.argv[1:]
flist.sort()

CH = 600
LSTBIN = .001
lsts = {}
raw_lst, raw_dat = [], []
for f in flist:
    print f
    uv = a.miriad.UV(f)
    uv.select('antennae', 0, 0)
    #uv.select('antennae', 1, 5)
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        #d /= gain(t)
        lst = n.round(uv['lst'] / LSTBIN) * LSTBIN
        try: lsts[lst].append(d[CH])
        except(KeyError): lsts[lst] = [d[CH]]
        raw_lst.append(uv['lst'])
        raw_dat.append(d[CH])
raw_lst = n.array(raw_lst)
raw_dat = n.array(raw_dat)
lst_t = lsts.keys()
lst_t.sort()
lst_t = n.array(lst_t)
for t in lst_t: lsts[t] = n.array(lsts[t])
dat_med = n.array([n.median(lsts[t]) for t in lst_t])
dat_mea = n.array([n.mean(lsts[t]) for t in lst_t])
dat_sig = n.array([n.sqrt(n.median(n.abs(lsts[t] - dat_med[i])**2)) \
        for i,t in enumerate(lst_t)])

# Rebin for smoother sigma
SIGBIN = .01
nbins = int(2*n.pi/SIGBIN)
swgts,sbins = n.histogram(lst_t, range=(0,2*n.pi), bins=nbins)
svals,sbins = n.histogram(lst_t, weights=dat_sig, range=(0,2*n.pi), bins=nbins)
sig_bin = svals/swgts.clip(1,n.inf)

# Do weighted average of data
dat_wgt_mea = []
for i,t in enumerate(lst_t):
    ind = n.floor(t / SIGBIN).clip(0,sig_bin.size-1).astype(n.int)
    wgts = n.exp(-2 * n.abs(lsts[t] - dat_med[i])**2 / sig_bin[ind])
    dat_wgt_mea.append(n.sum(lsts[t] * wgts) / n.sum(wgts))

p.subplot(211)
p.plot(lst_t, n.abs(dat_med), 'k.')
p.plot(lst_t, n.abs(dat_mea), 'r.')
p.plot(lst_t, n.abs(dat_wgt_mea), 'b')
#p.semilogy(raw_lst, n.abs(raw_dat), 'g.')
p.ylim(2,6)
#p.subplot(412)
#p.plot(lst_t, n.real(dat_med), 'k.')
#p.plot(lst_t, n.real(dat_mea), 'r.')
#p.plot(lst_t, n.real(dat_wgt_mea), 'b')
##p.plot(raw_lst, n.real(raw_dat), 'g.')
##p.ylim(-.25,.25)
#p.subplot(413)
#p.plot(lst_t, n.imag(dat_med), 'k.')
#p.plot(lst_t, n.imag(dat_mea), 'r.')
#p.plot(lst_t, n.imag(dat_wgt_mea), 'b')
##p.plot(raw_lst, n.imag(raw_dat), 'g.')
##p.ylim(-.25,.25)
p.subplot(212)
p.semilogy(lst_t, dat_sig, 'r.')
p.semilogy(.5*(sbins[1:] + sbins[:-1]), sig_bin, 'k')
#p.xlim(0,.15)
p.ylim(.0001,1)
p.show()
