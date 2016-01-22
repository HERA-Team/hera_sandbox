#! /usr/bin/env python
import aipy as a, numpy as n, sys, optparse, glob, ipdb, ephem, capo
from pylab import *
o=optparse.OptionParser()
o.set_usage("plot_lst_var.py [options]")
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True)
opts,args=o.parse_args(sys.argv[1:])
rmbls=[]
def get_data(filenames, antstr, polstr, rmbls, verbose=False):
    # XXX could have this only pull channels of interest to save memory
    lsts, cnt, var, dat, flg = [], {}, {}, {}, {}
    if type(filenames) == 'str': filenames = [filenames]
    for filename in filenames:
        if verbose: print '   Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, antstr, polstr)
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            bl = a.miriad.ij2bl(i,j)
            if bl in rmbls: continue
            lst = uv['lst']
            if len(lsts) == 0 or lst != lsts[-1]: 
                lsts.append(lst)
                #var.append(uv['var'])
            if not dat.has_key(bl):
                 dat[bl],flg[bl],cnt[bl],var[bl] = [],[],[],[]
            dat[bl].append(d)
            flg[bl].append(f)
            cnt[bl].append(uv['cnt'])
            var[bl].append(uv['var'])
    return n.array(lsts),cnt,var, dat, flg
lsts,cnt,var,data,flgs = get_data(args,antstr=opts.ant,polstr=opts.pol,rmbls=rmbls)
for bl in var:
    var[bl] = n.ma.masked_where(var[bl]==0,var[bl])*20000
for bl in data.keys():
    figure(1)
    clf()
    suptitle(str(a.miriad.bl2ij(bl))+opts.pol)
    subplot(131)
    title('data')
    imshow(n.abs(data[bl]),aspect='auto')
    colorbar()
    subplot(132)
    title('variance')
    imshow(n.sqrt(var[bl]),aspect='auto')
    print n.sqrt(var[bl]).max()
    colorbar()
    subplot(133)
    title('data/variance')
    SNR = n.abs(data[bl])/n.sqrt(var[bl]) 
    imshow(n.log10(SNR),aspect='auto',vmax=n.log10(n.mean(SNR)*2))
    colorbar()
    show()

