#! /usr/bin/env python
import aipy as a, numpy as n, os, sys, glob
import optparse
o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])

aa = a.phs.ArrayLocation(('-30:43:17.5', '21:25:41.9')) # Karoo, ZAR, GPS. #elevation=1085m


for filename in sys.argv[1:]:
    print filename, '->', filename+'c'
    if os.path.exists(filename+'c'):
        print '    File exists... skipping.'
        continue
    file_jd = float('.'.join(filename.split('.')[1:-1]))
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'c', status='new')
    curtime = 0

    def mfunc(uv, p, d, f):
        global curtime
        crd,t,(i,j) = p
        p1,p2 = a.miriad.pol2str[uv['pol']]
        # prevent multiple entries arising from xy and yx on autocorrelations
        if i == j and (p1,p2) == ('y','x'): return p, None, None

        if i > j: i,j,d = j,i,n.conjugate(d)
        if i == 8: d = -d  # I think the dipole for this antenna is rotated 180 deg
        if j == 8: d = -d

        if t != curtime:
            aa.set_jultime(t)
            uvo['lst'] = uvo['ra'] = uvo['obsra'] = aa.sidereal_time()
            curtime = t

        return (crd,t,(i,j)),d,f

    override = {
        'lst': aa.sidereal_time(),
        'ra': aa.sidereal_time(),
        'obsra': aa.sidereal_time(),
        'latitud': aa.lat,
        'dec': aa.lat,
        'obsdec': aa.lat,
        'longitu': aa.long,
        'nchan': 1024,
        'nants': 64,
        'ngains': 128,
        'nspect0': 64,
        'telescop':'PAPER',
        'antpos': n.transpose(n.array([
            [0., 0., 0.], #0
            [0., 0., 0.], #1
            [0., 0., 0.], #2
            [0., 0., 0.], #3
            [0., 0., 0.], #4
            [0., 0., 0.], #5
            [0., 0., 0.], #6
            [0., 0., 0.], #7
            [0., 0., 0.], #8
            [0., 0., 0.], #9
            [0., 0., 0.], #10
            [0., 0., 0.], #11
            [0., 0., 0.], #12
            [0., 0., 0.], #13
            [0., 0., 0.], #14
            [0., 0., 0.], #15
            [0., 0., 0.], #16
            [0., 0., 0.], #17
            [0., 0., 0.], #18
            [0., 0., 0.], #19
            [0., 0., 0.], #20
            [0., 0., 0.], #21
            [0., 0., 0.], #22
            [0., 0., 0.], #23
            [0., 0., 0.], #24
            [0., 0., 0.], #25
            [0., 0., 0.], #26
            [0., 0., 0.], #27
            [0., 0., 0.], #28
            [0., 0., 0.], #29
            [0., 0., 0.], #30
            [0., 0., 0.], #31
            [0., 0., 0.], #32
            [0., 0., 0.], #33
            [0., 0., 0.], #34
            [0., 0., 0.], #35
            [0., 0., 0.], #36
            [0., 0., 0.], #37
            [0., 0., 0.], #38
            [0., 0., 0.], #39
            [0., 0., 0.], #40
            [0., 0., 0.], #41
            [0., 0., 0.], #42
            [0., 0., 0.], #43
            [0., 0., 0.], #44
            [0., 0., 0.], #45
            [0., 0., 0.], #46
            [0., 0., 0.], #47
            [0., 0., 0.], #48
            [0., 0., 0.], #49
            [0., 0., 0.], #50
            [0., 0., 0.], #51
            [0., 0., 0.], #52
            [0., 0., 0.], #53
            [0., 0., 0.], #54
            [0., 0., 0.], #55
            [0., 0., 0.], #56
            [0., 0., 0.], #57
            [0., 0., 0.], #58
            [0., 0., 0.], #59
            [0., 0., 0.], #60
            [0., 0., 0.], #61
            [0., 0., 0.], #62
            [0., 0., 0.], #63
        ])).flatten(),
    }

    uvo.init_from_uv(uvi, override=override)

    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='\nCORRECT: '+' '.join(sys.argv)+'\n')
    del(uvo)
