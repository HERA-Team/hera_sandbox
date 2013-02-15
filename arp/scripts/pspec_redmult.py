#! /usr/bin/env python
import aipy as a, numpy as n
import capo as C
import optparse, sys, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True)
o.add_option('-t', '--taps', type='int', default=1,
    help='Taps to use in the PFB.  Default 1, which instead uses windowed FFT')
opts,args = o.parse_args(sys.argv[1:])

NTAPS = opts.taps
if NTAPS > 1: PFB = True
else: PFB = False
WINDOW = 'blackman-harris'
#WINDOW = 'none'

# XXX Currently hardcoded for PSA898
A_ = [0,16,8,24,4,20,12,28]
B_ = [i+1 for i in A_]
C_ = [i+2 for i in A_]
D_ = [i+3 for i in A_]
ANTPOS = n.array([A_, B_, C_, D_])

def bl_index(bl):
    i,j = a.miriad.bl2ij(bl)
    return i * 32 + j

# Get a dict of all separations and the bls that contribute
bl2sep = {}
sep2bl = {}
for ri in range(ANTPOS.shape[0]):
    for ci in range(ANTPOS.shape[1]):
        for rj in range(ANTPOS.shape[0]):
            for cj in range(ci,ANTPOS.shape[1]):
                if ri >= rj and ci == cj: continue # exclude repeat +/- listings of certain bls
                sep = a.miriad.ij2bl(rj-ri, cj-ci)
                i,j = ANTPOS[ri,ci], ANTPOS[rj,cj]
                bl = a.miriad.ij2bl(i,j)
                if i > j: i,j,sep = j,i,-sep
                bl2sep[bl] = sep
                sep = n.abs(sep)
                sep2bl[sep] = sep2bl.get(sep,[]) + [bl]

uv = a.miriad.UV(args[0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
sdf = uv['sdf']
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
del(uv)

afreqs = freqs.take(chans)

fq = n.average(afreqs)
z = C.pspec.f2z(fq)
if PFB:
    # XXX unsure how much of a BW modification a windowed PFB needs.  I think not much...
    B = sdf * afreqs.size / NTAPS
else:
    B = sdf * afreqs.size / C.pfb.NOISE_EQUIV_BW[WINDOW]
bm = n.polyval(C.pspec.DEFAULT_BEAM_POLY, fq)
scalar = C.pspec.X2Y(z) * bm * B
print fq, z, B

#cen_fqs = n.arange(.115,.190,.005)
#cen_fqs = n.array([.150])
#kwargs = {'cen_fqs':cen_fqs,'B':B, 'ntaps':NTAPS, 'window':WINDOW, 'bm_fqs':afreqs.clip(.120,.190)}
#window = a.dsp.gen_window(freqs.size, window=WINDOW)

for filename in args:
    outfile = filename + '.pspec'
    print filename,'->',outfile
    if os.path.exists(outfile):
        print '    File exists.  Skipping...'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi, override={'nchan':afreqs.size, 'sfreq':afreqs[0]})
    uvo._wrhd('history', uvi['history'] + 'PSPEC: ' + ' '.join(sys.argv) + '\n')
    #uvo.add_var('k3pk_fq', 'r')
    #uvo.add_var('k3pk_wgt', 'r')
    
    a.scripting.uv_selector(uvi, opts.ant, opts.pol)

    _Tlist,_Wlist,curtime = {},{},None
    for (crd,t,(i,j)),d,f in uvi.all(raw=True):
        if t != curtime:
            #print t
            uvo.copyvr(uvi)
            for sep,bls in sep2bl.items():
                for cnt,bl0 in enumerate(bls):
                    if not _Tlist.has_key(bl0): continue
                    for bl1 in bls[cnt:]: # this includes "auto-pspecs"
                        if not _Tlist.has_key(bl1): continue
                        pk = scalar * _Tlist[bl0] * n.conj(_Tlist[bl1])
                        uvo.write((crd,curtime,(bl_index(bl0),bl_index(bl1))), pk, n.zeros(pk.shape, dtype=n.int))
            # Clear the current pspec data and start a new integration
            _Tlist,_Wlist = {},{}
            curtime = t

        bl = a.miriad.ij2bl(i,j)
        sep = bl2sep[bl]
        if sep < 0: d,sep = n.conj(d),-sep

        d,f = d.take(chans), f.take(chans)
        w = n.logical_not(f).astype(n.float)
        Trms = d * C.pspec.jy2T(afreqs)
        if PFB:
            _Trms = C.pfb.pfb(Trms, window=WINDOW, taps=NTAPS, fft=n.fft.ifft)
            _Wrms = C.pfb.pfb(w   , window=WINDOW, taps=NTAPS, fft=n.fft.ifft)
        else:
            window = a.dsp.gen_window(Trms.size, WINDOW)
            _Trms = n.fft.ifft(window * Trms)
            _Wrms = n.fft.ifft(window * w)
        gain = n.abs(_Wrms[0])
        if gain > 0:
            _Tcln, info = a.deconv.clean(_Trms, _Wrms, tol=1e-9, maxiter=100, stop_if_div=False, verbose=False)
            #_Tcln, info = a.deconv.clean(_Trms, _Wrms, tol=1e-9)
            _Trms = _Tcln + info['res'] / gain
        _Trms = n.fft.fftshift(_Trms)
        _Wrms = n.fft.fftshift(_Wrms)

        _Tlist[bl] = _Trms
        _Wlist[bl] = _Wrms

    # Gotta do this one last time to catch the last integration.
    for sep,bls in sep2bl.items():
        for cnt,bl0 in enumerate(bls):
            if not _Tlist.has_key(bl0): continue
            for bl1 in bls[cnt:]: # this includes "auto-pspecs"
                if not _Tlist.has_key(bl1): continue
                pk = scalar * _Tlist[bl0] * n.conj(_Tlist[bl1])
                uvo.write((crd,curtime,(bl_index(bl0),bl_index(bl1))), pk, n.zeros(pk.shape, dtype=n.int))
    del(uvi); del(uvo)

