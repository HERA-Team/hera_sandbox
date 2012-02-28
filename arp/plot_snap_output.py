#! /usr/bin/env python
import numpy as n, pylab as p, sys, optparse, aipy as a

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, src=True)
opts, args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal, .1, .150, 1)
#RMANTS = [11]
RMANTS = []

symb = '-'
symb2 = ':'
state = {'time':None, 'object':None}
data = {'bms':{}} ; NPLTS = 3
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
        if True:
            aa.set_params(prms)
            cat.set_params(prms)
            aa.set_jultime(times[-1])
            lsts.append(aa.sidereal_time())
            cat.compute(aa)
            for srcname in cat:
              src = cat[srcname]
              if not data['bms'].has_key(srcname): data['bms'][srcname] = {}
              for o in prms:
                try: _o = int(o)
                except(ValueError): continue
                g = aa[_o].passband()
                b = aa[_o].bm_response(src.get_crds('eq',ncrd=3))
                if not data['bms'][srcname].has_key(o):
                    data['bms'][srcname][o] = []
                data['bms'][srcname][o].append(src.jys * n.abs(g*b)**2)
        else:
            aa.set_jultime(times[-1])
            lsts.append(aa.sidereal_time())
            cat.compute(aa)
            for srcname in cat:
              src = cat[srcname]
              if not data['bms'].has_key(srcname): data['bms'][srcname] = []
              b = aa[0].bm_response(src.get_crds('top',ncrd=3))**2
              data['bms'][srcname].append(b)  
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

lsts = n.array(lsts)
lsts = n.where(lsts < lsts[0], lsts + 2*n.pi, lsts)
dlst = n.concatenate([[0], lsts[1:] - lsts[:-1]])
mask = n.where(dlst < 0, 1, 0)
lsts = n.ma.array(lsts, mask=mask)
times = n.array(times)

for o in data:
    if o == 'bms':
        for srcname in data[o]:
          if type(data[o][srcname]) == dict:
            for ant in data[o][srcname]:
              data[o][srcname][ant] = n.array(data[o][srcname][ant]).squeeze()
              data[o][srcname][ant] = data[o][srcname][ant][:times.size]
          else:
              data[o][srcname] = n.array(data[o][srcname]).squeeze()
    else:
      for prm in data[o]:
        data[o][prm] = n.array(data[o][prm]).squeeze()
        data[o][prm] = data[o][prm][:times.size]

avg_ant, NANT = 0, 0
for o in data.keys():
    try: _o = int(o)
    except: continue
    avg_ant += data[o]['amp']
    NANT += 1
if NANT > 0:
    avg_ant /= NANT
    #if True: avg_ant = data['0']['amp']
    if True: avg_ant = data['2']['amp']

for o in data:
    try: _o = int(o)
    except(ValueError): _o = None
    if o == 'bms':
        p.subplot(NPLTS, 1, 3)
        for srcname in data[o]:
          if type(data[o][srcname]) == dict:
            for ant in data[o][srcname]:
              if int(ant) in RMANTS: continue
              p.plot(lsts, data[o][srcname][ant], symb2, label=str(ant))
          else:
              #p.plot(lsts, data[o][srcname], symb2, label=srcname)
              p.plot(lsts, data[o][srcname]/data[o]['cyg'], symb2, label=srcname)
    elif _o != None:
        p.subplot(NPLTS, 1, 2)
        if int(o) in RMANTS: continue
        #p.plot(data[o]['amp'], label=o)
        #p.plot(data[o]['amp']/avg_ant, label=o)
        d = data[o]['amp']/avg_ant
        d /= n.average(d[1:])
        if not o in ['0','2','3','4']: continue
        p.plot(lsts, d, symb, label=o)
    elif o in ['cas','cyg']:
        p.subplot(NPLTS, 1, 3)
        #p.plot(data[o]['jys'], label=o)
        #p.plot(lsts, data[o]['jys']/data['cyg']['jys'], symb, label=o)
        p.plot(lsts, data[o]['jys']*data['bms'][o]/data['cyg']['jys']/data['bms']['cyg'], symb, label=o)

p.subplot(NPLTS, 1, 1); p.plot(lsts, score, symb); p.ylim(.1,.3)
for i in range(NPLTS):
    p.subplot(NPLTS, 1, i+1)
    if False: p.legend()
    if False: p.xlim(5.2,6)
p.show()
