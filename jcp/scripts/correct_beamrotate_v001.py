#! /usr/bin/env python

'''A correct scripts for the beam rotation experiment.  Once again, this has been thrown together by JP so is likely terrible.'''

import aipy as a, numpy as n, os, optparse, sys, time

aa = a.phs.ArrayLocation(('38:25:59.24',  '-79:51:02.1')) # Green Bank, WV

antnumbers = { 
    0: {'x':1, 'y':13},
    1: {'x':1, 'y':13},
    2: {'x':3, 'y':17},
    3: {'x':3, 'y':17},
    4: {'x':6, 'y':18},
    5: {'x':6, 'y':18},
    6: {'x':8, 'y':19},
    7: {'x':8, 'y':19},
    8: {'x':9, 'y':27},
    9: {'x':9, 'y':27},
    10: {'x':10, 'y':15},
    11: {'x':10, 'y':28},
    12: {'x':12, 'y':29},
    13: {'x':12, 'y':29},
    14: {'x':28, 'y':31},
    15: {'x':16, 'y':31},
}

polnumbers = {
    0: {'x':'y', 'y':'y'},
    1: {'x':'x', 'y':'x'},
    2: {'x':'y', 'y':'y'},
    3: {'x':'x', 'y':'x'},
    4: {'x':'y', 'y':'y'},
    5: {'x':'x', 'y':'x'},
    6: {'x':'y', 'y':'y'},
    7: {'x':'x', 'y':'x'},
    8: {'x':'y', 'y':'y'},
    9: {'x':'x', 'y':'x'},
    10: {'x':'y', 'y':'x'},
    11: {'x':'x', 'y':'x'},
    12: {'x':'y', 'y':'y'},
    13: {'x':'x', 'y':'x'},
    14: {'x':'y', 'y':'y'},
    15: {'x':'x', 'y':'x'},
}

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])
orbchan = n.arange(278,287)


for filename in args:
    print filename, '->', filename+'c'
    if os.path.exists(filename+'c'):
        print '    File exists... skipping.'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'c', status='new')
    curtime = 0
    dat = []
    def mfunc(uv, p, d):
        global curtime, dat
        crd,t,(i,j) = p
        t += -2.0/24 #timezone error?
        p1,p2 = a.miriad.pol2str[uv['pol']]
        #print i, j, p1,p2

        oldi=i
        oldj=j
        oldp1=p1
        oldp2=p2
        
        if i == j and (p1,p2) == ('y','x'):
            return p, None

        if t!= curtime:
            aa.set_jultime(t)
            uvo['lst'] = aa.sidereal_time()
            uvo['ra'] = aa.sidereal_time()
            uvo['obsra'] = aa.sidereal_time()
            curtime = t
        #uvo['lst']=float(cntr)
        #print cntr,uvo['lst']
        if (j-i)>8:
            d = n.conjugate(d);
            if p1 != p2:
                p1,p2 = p2,p1
                d = n.conjugate(d)
        if (p1,p2)==('y','x'):
            d = n.conjugate(d)

        tp1 = polnumbers[i][p1]
        tp2 = polnumbers[j][p2]
        temp_pol = tp1 + tp2
        
        i = antnumbers[i][p1]
        j = antnumbers[j][p2]

        #if temp_pol == 'xx' and i==1 and j==1:
            #print oldi,oldj,oldp1,oldp2

        uvo['pol'] = a.miriad.str2pol[temp_pol]
        #print a.miriad.str2pol[temp_pol],uvo['pol']
        #print temp_pol, a.miriad.pol2str[uvo['pol']]
        #if a.miriad.str2pol[temp_pol] != uvo['pol']: print 'awooooogah!'
        #print '--------------'
        
        if i > j:
            i,j = j,i
            d = n.conjugate(d)

        p = crd,t,(i,j)
        #cntr +=1
        if temp_pol == 'xx' and i==1 and j==1:
            print "DID NOT blank", i,j,temp_pol
            return p,d
        elif temp_pol in ['yy','xx','xy','yx']:
            return p,d
        else:
            return None, None
            
    override = {
        'lst': aa.sidereal_time(),
        'ra': aa.sidereal_time(),
        'obsra': aa.sidereal_time(),
        'sdf': 97.65e-6,
        'latitud': aa.lat,
        'dec': aa.lat,
        'obsdec': aa.lat,
        'longitu': aa.long,
        'sfreq':0.10009765625,
        'freq':0.10009765625,
        'inttime':200000000./1024/128/4096,
        'nants': 33,
        'ngains': 64,
        'nspect0': 32,
        'pol':a.miriad.str2pol['xy'],
        'telescop':'PAPER',
        'bandpass': n.resize(uvi['bandpass'], (33,1024)).flatten(),
        'antpos': n.transpose(n.array([
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.],
        ])).flatten(),
    }

    uvo.init_from_uv(uvi, override=override)
    uvo.pipe(uvi, mfunc=mfunc,
        append2hist='Rewrote sdf/sfreq/freq/inttime to be correct\nSet initial antenna positions to ZERO.\nChanged times by -2 hours.\nFixed? correlator conjugation errors.')
    #n.save('xx11',dat)
    del(uvo)
