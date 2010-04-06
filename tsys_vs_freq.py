#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, re, glob, sys

N_ant = 8
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
        60 : 3170.0, 65 : 3351.4, 70 : 3361.2, 75 : 2698.6,
        80 : 3178.6, 85 : 2082.2, 90 : 2816.0, 95 : 2132.6,
        100 : 2993.8, 105 :  333.8, 110 : 2503.6, 115 : 3062.8,
        120 : 2642.4, 125 :   54.6, 130 : 3309.2, 135 : 1655.8,
        140 :    0.0, 145 : 3036.0, 150 : 3390.2, 155 : 3182.8,
        160 : 3426.4, 165 : 3214.4, 170 : 2634.2, 175 : 3317.4,
    },
    '11:20_40:00' : {
        60 : 3565.5, 65 : 3727.2, 70 : 3755.8, 75 : 3115.5,
        80 : 3529.0, 85 : 2270.7, 90 : 3199.7, 95 : 2352.5,
        100 : 2937.0, 105 :  697.0, 110 : 2933.2, 115 : 3453.8,
        120 : 2456.0, 125 :  615.3, 130 : 3696.3, 135 : 1538.5,
        140 :  288.2, 145 : 3453.7, 150 : 3707.3, 155 : 3583.8,
        160 : 3649.8, 165 : 3597.0, 170 : 3076.3, 175 : 3059.5,
    },
    '12:00_30:00' : {
        60 : 3246.3, 65 : 3373.0, 70 : 3388.3, 75 : 2838.8,
        80 : 3242.5, 85 : 2264.7, 90 : 2915.0, 95 : 1973.7,
        100 : 2653.2, 105 :  810.3, 110 : 2853.2, 115 : 3161.2,
        120 : 2236.2, 125 :  619.8, 130 : 3350.0, 135 : 1309.0,
        140 :  179.8, 145 : 3035.0, 150 : 3363.8, 155 : 3268.0,
        160 : 3280.5, 165 : 3241.0, 170 : 2784.5, 175 : 2712.2,
    },
    '12:00_40:00' : {
        60 : 3621.3, 65 : 3758.0, 70 : 3775.0, 75 : 3142.3,
        80 : 3605.0, 85 : 2533.5, 90 : 3267.0, 95 : 2191.0,
        100 : 2943.3, 105 :  886.5, 110 : 3145.2, 115 : 3515.3,
        120 : 2494.7, 125 :  701.7, 130 : 3726.0, 135 : 1466.7,
        140 :  210.0, 145 : 3387.0, 150 : 3743.0, 155 : 3640.2,
        160 : 3656.0, 165 : 3614.3, 170 : 3101.3, 175 : 3006.3,
    },
    #'11:20_30:00': {
    #    60 : 1949.8, 65 : 2022.0, 70 : 2019.8, 75 : 1561.4,
    #    80 : 1903.4, 85 : 1239.6, 90 : 1681.2, 95 : 1216.6,
    #    100 : 1815.8, 105 : 100.8, 110 : 1469.2, 115 : 1857.8,
    #    120 : 1556.2, 125 : 47.6, 130 : 1972.2, 135 : 1055.4,
    #    140 : 0.0, 145 : 1889.0, 150 : 2020.0, 155 : 1947.0,
    #    160 : 2031.4, 165 : 1907.0, 170 : 1601.4, 175 : 2009.8,
    #},
    #'11:20_40:00' : {
    #    60 : 2268.6, 65 : 2344.6, 70 : 2344.6, 75 : 1804.4,
    #    80 : 2212.2, 85 : 1413.8, 90 : 1932.4, 95 : 1345.4,
    #    100 : 2052.2, 105 : 136.8, 110 : 1716.6, 115 : 2161.2,
    #    120 : 1789.2, 125 : 62.4, 130 : 2285.0, 135 : 1214.4,
    #    140 : 0.0, 145 : 2175.6, 150 : 2343.2, 155 : 2259.2,
    #    160 : 2355.6, 165 : 2212.0, 170 : 1856.0, 175 : 2310.6,
    #},
    #'12:00_30:00' : {
    #    60 :  2052.6, 65 :  2111.0, 70 :  2107.4, 75 :  1616.2,
    #    80 :  2031.4, 85 :  1439.0, 90 :  1770.4, 95 :  1115.2,
    #    100 :  1795.8, 105 :   177.6, 110 :  1724.2, 115 :  1980.0,
    #    120 :  1595.6, 125 :    84.6, 130 :  2067.6, 135 :  1019.4,
    #    140 :     0.0, 145 :  1891.2, 150 :  2111.0, 155 :  2049.8,
    #    160 :  2117.0, 165 :  1966.2, 170 :  1669.8, 175 :  2018.8,
    #},
    #'12:00_40:00': {
    #    60 : 2352.6, 65 : 2417.8, 70 : 2414.4, 75 : 1840.4,
    #    80 : 2317.4, 85 : 1656.4, 90 : 2048.6, 95 : 1269.4,
    #    100 : 2039.2, 105 : 201.8, 110 : 1952.0, 115 : 2259.4,
    #    120 : 1820.8, 125 : 98.0, 130 : 2362.2, 135 : 1170.6,
    #    140 : 0.0, 145 : 2163.4, 150 : 2415.6, 155 : 2348.4,
    #    160 : 2424.4, 165 : 2253.2, 170 : 1912.8, 175 : 2293.2,
    #},
}

def ch2freq(chan, sfreq=0.1212890625, sdf=0.00029296875):
    return (sfreq + chan*sdf) * 1e9

def T(Jy, chan=162.5, sfreq=0.1212890625, sdf=0.00029296875):
    freq = ch2freq(chan, sfreq=sfreq, sdf=sdf)
    lam = a.const.c / freq
    bm = 10**n.polyval(bm_poly, n.log10(freq))
    return Jy * 1e-23 / bm * lam**2 / (2 * a.const.k)

freqs, tsys, tsky, tnos = {}, {}, {}, {}
for src in t_3day:
    freqs[src], tsys[src], tsky[src], tnos[src] = [], [], [], []

for pfx in prefixes:
    f2a = pfx+'00a.bim.fits'
    f2b = pfx+'00b.bim.fits'
    ch1, ch2 = map(float, f_re.match(pfx).groups())
    ch = (ch1 + ch2) / 2.
    freq = ch2freq(ch)
    pb = n.polyval(pb_poly, freq*1e-9)
    bm = 10**n.polyval(bm_poly, n.log10(freq))

    srms0,srms1 = [], []
    d2a, kwds = a.img.from_fits(f2a); d2a = d2a.squeeze()
    d2b, kwds = a.img.from_fits(f2b); d2b = d2b.squeeze()
    src = kwds['object']
    if t_3day[src][ch1] < 1000: continue
    t = t_3day[src][ch1] * T_int
    d1 = (d2a + d2b) / 2.
    d2 = (d2a - d2b) / 2.

    #px_rms1 = n.sqrt(n.average(d1[466:533,466:533]**2))
    #px_rms2 = n.sqrt(n.average(d2[466:533,466:533]**2))
    px_rms1 = n.std(d1[466:533,466:533])
    px_rms2 = n.std(d2[466:533,466:533])

    print '-----------------------------------------------------'
    print pfx, src
    print 'FREQ =', freq/1e6
    #print 'SKY RMS1:', px_rms1
    #print 'NOISE RMS2:', px_rms2
    #print 'T =', T(px_rms1, chan=ch)
    #print 'dT =', T(px_rms2, chan=ch)
    #print 't_int =', t
    #print 'PB =', pb
    #print 'BM =', bm
    Tsys = T(px_rms2, chan=ch) * n.sqrt(BW  * t * N_ant * (N_ant-1)) * bm / pb
    #print 'Tsys:', Tsys

    freqs[src].append(freq)
    tsys[src].append(Tsys)
    tsky[src].append(T(px_rms1, chan=ch))
    tnos[src].append(T(px_rms2, chan=ch))

tsys_tot, tsky_tot, tnos_tot = 0, 0, 0
fqs = None
for src in t_3day:
    f = n.array(freqs[src])
    print src
    o = n.argsort(f)
    f = f.take(o)
    if fqs is None: fqs = f
    print f
    print fqs
    print '-'*70
    assert(n.all(fqs == f))
    tsys_f = n.array(tsys[src]).take(o)
    tsky_f = n.array(tsky[src]).take(o)
    tnos_f = n.array(tnos[src]).take(o)
    tsys_tot += tsys_f
    tsky_tot += tsky_f
    tnos_tot += tnos_f**2
    #tbalun_vs_fq = n.polyval(tbalun_poly, freqs)
    #tsky_vs_fq = n.sqrt(tsys_vs_fq**2 - tbalun_vs_fq**2)
    #p.plot(f/1e6, tsys_f, '.-', label='Tsys')
    #p.errorbar(f/1e6, tsky_f, tnos_f, fmt='-.', label='Tsky')
tsys_tot /= len(t_3day)
tsky_tot /= len(t_3day)
tnos_tot /= len(t_3day)
tnos_tot = n.sqrt(tnos_tot)

def tsky_sync(a, b): return a * (n.arange(135.,175)/ 150.)**b
#p.plot(n.arange(135.,175), tsky_sync(240,-2.5), '-', label='Tsync')
##ax = p.gca(); ax.set_xscale('log'); ax.set_yscale('log')
#p.xlim(135, 175)
#p.ylim(0, 350)
##p.legend()
#p.show()

p.plot(n.arange(135.,175), tsky_sync(240,-2.5), 'k-.', label='Tsync')
#p.plot(fqs/1e6, tsys_tot, 'k.-', label='Tsys')
p.errorbar(fqs/1e6, tsys_tot, 2/66**.5*tsys_tot, fmt='k.-', label='Tsys')
p.errorbar(fqs/1e6, tsky_tot, 2*tnos_tot, fmt='k:.', label='Tsky')
#p.plot(fqs/1e6, tnos_tot, 'k^-', label='Tpx')
p.errorbar(fqs/1e6, tnos_tot, 2/66**.5*tnos_tot, fmt='k^-', label='Tpx')
ax = p.gca(); ax.set_xscale('log'); ax.set_yscale('log')
p.xlim(135, 175)
p.ylim(1e0, 1e3)
p.xticks(n.arange(135,180,5), n.arange(135,180,5))
p.xlabel('Frequency (MHz)')
p.ylabel('Temperature (K)')
p.show()
