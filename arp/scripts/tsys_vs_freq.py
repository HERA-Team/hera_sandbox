#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, re, glob, sys

N_ant = 8
#N_ant = 12
BW = 1.46e6 # = 5 channels
T_int = 14.32
bm_poly = [ -2.02131219,  11.89480783]
pb_poly = [1.02854332e+09, -9.49707493e+08, 3.64775002e+08, 
    -7.46038156e+07, 8.56951433e+06, -5.24246222e+05, 1.33464786e+04]

tbalun_fq = [72.4, 76.3, 80.2, 84.1, 88, 91.9, 95.8, 99.7, 103.6,
    107.5, 111.4, 115.3, 119.2, 123.1, 127, 130.9, 134.8, 138.7, 142.6,
    146.5, 150.4, 154.3, 158.2, 162.1, 166, 169.9, 173.8, 177.7, 181.6,
    185.5, 189.4, 193.3, 197.2, 201.1, 205, 208.9, 212.8, 216.7, 220.6, 224.5,]
tbalun_ms = [112.31, 119.67, 115.18, 118.22, 120.56, 123.51, 127.33,
    125.86, 134.3, 125.58, 122.95, 121.03, 121.22, 118.25, 119.32, 114.43,
    112.62, 112.02, 109.84, 108.82, 108.74, 108.84, 108.1, 109, 109.5, 110.4,
    110.6, 112.16, 117.13, 118.49, 121.21, 122.62, 126.24, 131.38, 131.65, 
    134.88, 136.84, 138.07, 143.19, 141.74,]
tbalun_poly = n.polyfit(n.array(tbalun_fq)*1e6, n.array(tbalun_ms), deg=6)

f_re = re.compile(r'.*_c(\d+)_(\d+)_.*')

prefixes = {}
for filename in sys.argv[1:]:
    pfx = filename[:filename.find('_0')+1]
    print filename, pfx
    prefixes[pfx] = None

t_3day = {
    '11:20_30:00' : {
    60 :     1646.8, 65 :     1711.5, 70 :     1675.6, 75 :     1339.4,
    80 :     1541.5, 85 :     1078.4, 90 :     1383.2, 95 :      948.8,
    100 :     1357.6, 105 :      145.8, 110 :     1195.4, 115 :     1421.7,
    120 :     1238.6, 125 :       20.3, 130 :     1565.7, 135 :      787.7,
    140 :        0.0, 145 :     1464.5, 150 :     1641.2, 155 :     1538.1,
    160 :     1665.8, 165 :     1576.6, 170 :     1289.0, 175 :     1632.4,
    },
    '11:20_40:00' : {
    60 :     2146.3, 65 :     2239.0, 70 :     2203.1, 75 :     1766.7,
    80 :     2042.6, 85 :     1428.5, 90 :     1846.3, 95 :     1271.5,
    100 :     1827.7, 105 :      197.7, 110 :     1618.1, 115 :     1931.6,
    120 :     1686.2, 125 :       28.5, 130 :     2142.0, 135 :     1079.6,
    140 :        0.0, 145 :     2017.2, 150 :     2264.7, 155 :     2126.3,
    160 :     2307.9, 165 :     2188.6, 170 :     1790.4, 175 :     2267.7,
    },
    '12:00_30:00' : {
    60 :     1690.1, 65 :     1731.2, 70 :     1689.2, 75 :     1346.8,
    80 :     1576.1, 85 :     1118.8, 90 :     1358.0, 95 :      892.9,
    100 :     1333.8, 105 :      215.5, 110 :     1318.0, 115 :     1457.5,
    120 :     1251.2, 125 :       48.7, 130 :     1580.0, 135 :      740.3,
    140 :        0.0, 145 :     1419.9, 150 :     1652.7, 155 :     1579.3,
    160 :     1672.5, 165 :     1578.2, 170 :     1305.3, 175 :     1616.9,
    },
    '12:00_40:00' : {
    60 :     2200.1, 65 :     2264.1, 70 :     2220.4, 75 :     1776.6,
    80 :     2088.5, 85 :     1484.6, 90 :     1814.7, 95 :     1198.2,
    100 :     1796.6, 105 :      291.4, 110 :     1782.8, 115 :     1979.8,
    120 :     1703.7, 125 :       67.2, 130 :     2161.7, 135 :     1015.0,
    140 :        0.0, 145 :     1956.2, 150 :     2281.1, 155 :     2184.8,
    160 :     2317.7, 165 :     2191.3, 170 :     1814.3, 175 :     2248.4,
    },
}

def ch2freq(chan, sfreq=0.1212890625, sdf=0.00029296875):
    return (sfreq + chan*sdf) * 1e9

def T(Jy, bm, chan=162.5, sfreq=0.1212890625, sdf=0.00029296875):
    freq = ch2freq(chan, sfreq=sfreq, sdf=sdf)
    lam = a.const.c / freq
    return Jy * 1e-23 / bm * lam**2 / (2 * a.const.k)

def calc_bm_area(dbm, dra, ddec):
    if True:
        #print adbm.max(), n.average(adbm), dra * ddec
        return n.sqrt(n.sum(n.abs(dbm/dbm.max())**2 * dra * ddec))
        #return n.where(adbm > n.max(adbm)/2, dra*ddec, 0).sum()
        #return n.where(adbm >= 1./8, dra*ddec * adbm**2, 0).sum()
    else:
        _dbm = n.fft.fft2(dbm)
        _adbm = n.abs(_dbm)
        _adbm /= _adbm.max()
        #return n.where(_adbm >= 1./8, dra*ddec * _adbm, 0).sum()
        #return .4 / _adbm.sum()
        #return .4 / (28*250.)
        return .4 / n.where(_adbm >= 1./16, 1, 0).sum()

freqs, tsys, tsky, trms = {}, {}, {}, {}
for src in t_3day:
    freqs[src], tsys[src], tsky[src], trms[src] = [], [], [], []

for pfx in prefixes:
    print '='*70
    print pfx
    ch1, ch2 = map(float, f_re.match(pfx).groups())
    ch = (ch1 + ch2) / 2.
    freq = ch2freq(ch)
    pb = n.polyval(pb_poly, freq*1e-9)

    srms0,srms1 = [], []
    d2a, kwds = a.img.from_fits(pfx+'00a.dim.fits'); d2a = d2a.squeeze()
    d2b, kwds = a.img.from_fits(pfx+'00b.dim.fits'); d2b = d2b.squeeze()
    b2a, kwds = a.img.from_fits(pfx+'00a.dbm.fits'); b2a = b2a.squeeze()
    b2b, kwds = a.img.from_fits(pfx+'00b.dbm.fits'); b2b = b2b.squeeze()
    b1 = b2a + b2b
    d1 = d2a + d2b
    d2 = d2a - d2b
    bm_gain = b1.max()
    bm_area = calc_bm_area(b1, kwds['d_ra'] * a.img.deg2rad, kwds['d_dec'] * a.img.deg2rad)
    print bm_area, 10**n.polyval(bm_poly, n.log10(freq))
    src = kwds['object']
    #if t_3day[src][ch1] == 0: continue
    if t_3day[src][ch1] < 500: continue
    t = t_3day[src][ch1] * T_int
    px_rms1 = n.std(d1[466:533,466:533]) / bm_gain
    px_rms2 = n.std(d2[466:533,466:533]) / bm_gain

    print '-----------------------------------------------------'
    print pfx, src
    print 'FREQ =', freq/1e6
    print 'SKY RMS1:', px_rms1
    print 'NOISE RMS2:', px_rms2
    print 't_int =', t
    print 'AREA PB =', pb
    print 'GAIN BM =', bm_gain
    #print 'AREA BM =', bm_area
    Tsys = T(px_rms2, pb, chan=ch) * n.sqrt(2*BW*t * (N_ant*(N_ant-1)/2))
    print 'Tsys:', Tsys
    if False:
        _d1 = n.fft.fft2(d1)
        _d2 = n.fft.fft2(d2)
        _b1 = n.fft.fft2(b1)
        for cnt,dat in enumerate([_d1,_b1,_d1/_b1,_d2,_b1,_d2/_b1]):
            p.subplot(2,3,cnt+1)
            plot_dat = a.img.recenter(n.log10(n.abs(dat)), (500,500))
            if cnt in [2,5]: vmax,vmin = 2,0
            else: vmax,vmin = plot_dat.max(), plot_dat.max()-2
            p.imshow(plot_dat, origin='lower', vmax=vmax, vmin=vmin)
            p.colorbar(shrink=.3)
        p.show()
        import sys; sys.exit(0)
    elif True:
        _d1 = n.fft.fft2(d1)
        _d2 = n.fft.fft2(d2)
        _b1 = n.fft.fft2(b1)
        if False:
            n1_uv = n.sqrt((n.abs(_d1)**2).sum() / (n.abs(_b1)**2).sum())
            n2_uv = n.sqrt((n.abs(_d2)**2).sum() / (n.abs(_b1)**2).sum())
        else:
            n1_uv = n.abs(_d1).sum() / n.abs(_b1).sum()
            n2_uv = n.abs(_d2).sum() / n.abs(_b1).sum()
        Tsky = T(n1_uv, pb, chan=ch)
        Trms = T(n2_uv, pb, chan=ch)
    else:
        Tsky = T(px_rms1, 1., chan=ch) / bm_area
        Trms = T(px_rms2, 1., chan=ch) / bm_area
    print 'Tsky:', Tsky
    print 'Trms:', Trms

    freqs[src].append(freq)
    tsys[src].append(Tsys)
    tsky[src].append(Tsky)
    trms[src].append(Trms)

tsys_tot, tsky_tot, trms_tot = 0, 0, 0
wgt_tot = 0
fqs = None
for src in t_3day:
    f = n.array(freqs[src])
    print src
    o = n.argsort(f)
    f = f.take(o)
    if fqs is None: fqs = f
    print f/1e9
    print fqs/1e9
    print '-'*70
    assert(n.all(fqs == f))
    tsys_f = n.array(tsys[src]).take(o)
    tsky_f = n.array(tsky[src]).take(o)
    trms_f = n.array(trms[src]).take(o)
    wgt = 1./trms_f
    tsys_tot += tsys_f * wgt
    tsky_tot += tsky_f * wgt
    wgt_tot += wgt
    #trms_tot += trms_f**2
    trms_tot += trms_f * wgt
tsys_tot /= wgt_tot
tsky_tot /= wgt_tot
trms_tot /= wgt_tot
#trms_tot /= len(t_3day)

def tsky_sync(a, b): return a * (n.arange(135.,175)/ 150.)**b

p.plot(n.arange(135.,175), tsky_sync(240,-2.5), 'k-.', label='Tsync')
#p.plot(n.arange(135.,175), tsky_sync(240,-2.5)+n.polyval(tbalun_poly, n.arange(135.,175)*1e6), 'k-', label='Tsync+Tbalun')
p.errorbar(fqs/1e6, tsys_tot, 2/66**.5*tsys_tot, fmt='k.-', label='Tsys')
p.errorbar(fqs/1e6, tsky_tot, 2*trms_tot, fmt='k:.', label='Tsky')
p.errorbar(fqs/1e6, trms_tot, 2/66**.5*trms_tot, fmt='k^-', label='Tpx')
ax = p.gca(); ax.set_xscale('log'); ax.set_yscale('log')
p.xlim(135, 175)
p.ylim(1e-2, 5e3)
p.xticks(n.arange(135,180,5), n.arange(135,180,5))
p.xlabel('Frequency (MHz)')
p.ylabel('Temperature (K)')
p.show()
