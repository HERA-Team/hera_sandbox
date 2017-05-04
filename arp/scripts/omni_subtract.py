#! /usr/bin/env python

import omnical, aipy, numpy, capo
import pickle, optparse, os, sys

### Options ###
o = optparse.OptionParser()
o.set_usage('omni_apply.py [options] *uv')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o, pol=True, cal=True)
opts,args = o.parse_args(sys.argv[1:])

aa = aipy.cal.get_aa(opts.cal,numpy.array([.15]))
print 'Getting reds from calfile'
info = capo.omni.aa_to_info(aa, pols=['y'], ex_ants=[], crosspols=['yy'])
reds = info.get_reds()

def get_red(sep):
     for r in reds:
          if sep in r: return r

class Baselines:
    def __init__(self, data, aliases={}):
        self.data = data
        self.aliases = aliases
    def __getitem__(self, name):
        try: name = self.aliases[name]
        except(KeyError):
            try: name = self.aliases[name[::-1]][::-1]
            except(KeyError): pass
        try: return self.data[name]
        except(KeyError): return self.data[name[::-1]].conj()

### Read Data and Solutions ###
for filename in args:
    npzfile = '.'.join(filename.split('.')[:4]+['npz'])
    print 'Reading', filename, '(%s)' % npzfile
    newfile = filename + 'O'
    if os.path.exists(newfile):
        print '    %s exists.  Skipping...' % newfile
        continue
    m,g,v,_ = capo.omni.from_npz(npzfile)
    for pol in v:
        aliases = {}
        for bl in v[pol]:
            for blr in get_red(bl): aliases[blr] = bl
        v[pol] = Baselines(v[pol], aliases)
                
    rfi = capo.xrfi.omni_chisq_to_flags(m['chisq'])
    times = []
    
    def mfunc(uv,p,d,f): #loops over time and baseline
        global times #global list
        _,t,(a1,a2) = p
        p1,p2 = pol = aipy.miriad.pol2str[uv['pol']]
        if len(times) == 0 or times[-1] != t: times.append(t)
        ti = len(times) - 1 #time index
#        import IPython; IPython.embed()
        try:
            mdl = g[p1][a1][ti] * g[p2][a2][ti].conj() * v[pol][(a1,a2)][ti]
            return p, d-mdl, rfi[ti]
        except(KeyError):
            #print pol, a1, a2
            #if a1 != a2:
            #    import IPython; IPython.embed()
            return p, d, rfi[ti]
    
    uvi = aipy.miriad.UV(filename)
    uvo = aipy.miriad.UV(newfile,status='new')
    uvo.init_from_uv(uvi)
    print '    Saving', newfile
    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist='OMNICAL: ' + ' '.join(sys.argv) + '\n')
        
