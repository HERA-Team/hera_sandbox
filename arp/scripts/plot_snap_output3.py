#! /usr/bin/env python
import numpy as n, pylab as p, sys, optparse, aipy as a

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, src=True)
opts, args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal, .1, .150, 1)
# Not so good: 5,
ANTS = ['5','14']

symb = '-'
symb2 = ':'
state = {'time':None, 'object':None}
data = {'bms':{}}
prms = {}
times, lsts = [], []
score = []
srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
bm_resp = {}

for filename in args:
  print 'Reading', filename
  for L in open(filename).readlines():
    W = L.split()
    if len(W) == 0: continue
    elif W[0].startswith('Time'):
        state['time'] = float(W[1])
    elif W[0].startswith('Score'):
        sys.stdout.write('.'); sys.stdout.flush()
        try:
            scr = float(W[1])
            if scr == 0: continue
            score.append(scr)
        except(IndexError): continue
        times.append(state['time'])
        #print prms
        aa.set_params(prms)
        cat.set_params(prms)
        aa.set_jultime(times[-1])
        lsts.append(aa.sidereal_time())
        cat.compute(aa)
        for srcname in cat:
          src = cat[srcname]
          if not data['bms'].has_key(srcname): data['bms'][srcname] = {}
          for o in prms.keys():
            if not o in ANTS + ['0']: continue
            _o = int(o)
            g = aa[_o].passband()
            b = aa[_o].bm_response(src.get_crds('top',ncrd=3))
            if not data['bms'][srcname].has_key(o):
                data['bms'][srcname][o] = []
            #if srcname == 'cyg': print _o, b, aa[_o].amp, aa[_o].beam.cyg_cas_ratio
            data['bms'][srcname][o].append(src.jys * n.abs(g*b)**2)
        state['score'] = float(W[1])
        state['time'] = None
        state['object'] = None
        prms = {}
    elif L.startswith(' ') and len(W) == 1 and state['time'] != None:
        state['object'] = W[0]
    elif L.startswith('   ') and state['object'] != None:
        o,prm = state['object'], W[0]
        if not prms.has_key(o): prms[o] = {}
        if not data.has_key(o): data[o] = {}
        if not data[o].has_key(prm): data[o][prm] = []
        prms[o][prm] = map(float, W[1:])
        data[o][prm].append(map(float, W[1:]))
        if len(prms[o][prm]) == 1: prms[o][prm] = prms[o][prm][0]
  print 

lsts = n.array(lsts)
lsts = n.where(lsts < lsts[0], lsts + 2*n.pi, lsts)
dlst = n.concatenate([[0], lsts[1:] - lsts[:-1]])
mask = n.where(dlst < 0, 1, 0)
lsts = n.ma.array(lsts, mask=mask)
times = n.array(times)

for o in data:
    if o == 'bms':
        for srcname in data[o]:
            for ant in data[o][srcname]:
                data[o][srcname][ant] = n.array(data[o][srcname][ant]).squeeze()
                data[o][srcname][ant] = data[o][srcname][ant][:times.size]
    else:
      for prm in data[o]:
        data[o][prm] = n.array(data[o][prm]).squeeze()
        data[o][prm] = data[o][prm][:times.size]

_dat = data['bms']
for srcname in _dat:
    if srcname == 'cyg': off=1
    else: off=0
    for ant in _dat[srcname]:
        if not ant in ANTS: continue
        p.subplot(3, 1, 1)
        p.plot(lsts, _dat[srcname][ant]/_dat['cyg'][ant], symb, label=srcname+ant)
        p.subplot(3, 1, 2)
        #p.plot(lsts, _dat[srcname][ant]/_dat[srcname]['0'], symb, label=srcname+ant)
        #p.plot(lsts, data[ant]['amp'], symb, label='amp'+ant)
        #p.plot(lsts, data[ant]['amp']*data[ant]['cyg_cas_ratio'], symb, label='rat'+ant)
        p.plot(lsts, _dat[srcname][ant], symb, label=srcname+ant)

p.subplot(3, 1, 3); p.plot(lsts, score, symb); p.ylim(.1,.3)
for i in range(3):
    p.subplot(3, 1, i+1)
    if False: p.legend()
    if False: p.xlim(5.2,6)
p.show()
