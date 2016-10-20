#!/usr/bin/env python
"""
Plot phase ratios of redundant baselines, optionally grouped by f-engine.
"""
import capo as C, numpy as np, matplotlib.pyplot as plt
import sys, subprocess, optparse, aipy

o = optparse.OptionParser()
o.set_usage('plot_phase_ratios.py [options] *.uvcRRE')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o, cal=True, pol=True)
o.add_option('-f', '--fengine', dest='fno', default=None,\
             help='Comma-separated f-engine numbers to source baselines from.')
o.add_option('-m', '--mfile', dest='mfile', default='/home/saulkohn/ReposForCanopy/capo/sak/data/fx_to_ant.dat',\
             help="Map file for fengine-to-antenna. Assumed mapping ['fengine', 'plate', 'rx', 'stations', 'ant'] with '->' separators")             
o.add_option('-s', '--seps', dest='seps', default="0,2",\
             help='Separation type to source baselines from. e.g. "0,1;0,2"')
o.add_option('-r','--ref', dest='ref', default=None,\
             help='We can plot with reference to a single baseline if we want. Specify array index (int).')
o.add_option('-o', '--outpath', dest='outpath', default=None,\
             help='Outpath for images. If None, no images will be saved.')
o.add_option('-v', '--verbose', dest='verb', action='store_true', help='Toggle verbosity')
opts,args = o.parse_args(sys.argv[1:])

if opts.verb:
    for i in range(1,9):
        ants_str = C.sak.get_ants(opts.mfile,pol=opts.pol,search='f%i'%i)
        print 'fengine %i: %s'%(i,ants_str) 

def bl2tup(bl):
    a1,a2 = map(int,bl.split('_'))
    return (a1,a2)

# prep output folder
if not opts.outpath is None:
    if not os.path.exists(opts.outpath):
        os.mkdir(opts.outpath)

# get info object from cal file
aa = aipy.cal.get_aa(opts.cal,np.array([.150]))
info = C.omni.aa_to_info(aa)

# get list of redundant baselines (strings)
proc = subprocess.Popen('grid2ant.py -C %s --seps="%s"'%(opts.cal,opts.seps), shell=True, stdout=subprocess.PIPE)
bls_str = proc.communicate()[0].rstrip()
red_bls_arr = bls_str.split(',') #['0_1','1_2',...]
"""
# parse f-engine (ints)
if not opts.fno is None: fengines = map(int,opts.fno.split(','))
else: fengines = range(1,9)

N_red_bls_per_fengine = []
for cnt,i in enumerate(fengines):
    if opts.verb: print 'fengine %i'%i
    ants_str = C.sak.get_ants(opts.mfile,pol=opts.pol,search='f%i'%i)
    ants = map(int,ants_str.split(',')) #[0,1,2,...]
    f_bls = []
    for j in ants:
        for k in ants:
            if j!=k:
                f_bls.append('%i_%i'%(j,k))
    pltbls = []
    for bl in f_bls:
        if bl in red_bls_arr: 
            pltbls.append(bl)
    
    N_red_bls_per_fengine.append(len(pltbls))
    
    pltbls_str=','.join(pltbls)
    if opts.verb: print pltbls_str
    #actually get data
    tinfo,data,flags = C.arp.get_dict_of_uv_data(args,pltbls_str,opts.pol)
    odata = C.sak.order_data(data,info)
    if not opts.ref is None: C.sak.plot_phase_ratios(odata, ref=int(opts.ref))
    else: C.sak.plot_phase_ratios(odata)
    if not opts.outpath is None: plt.savefig(opts.outpath+'/fengine_%i_refbl=%s.png'%(i,pltbls_str[int(opts.ref)]))
#plt.show()
plt.close()
"""
nbls = len(red_bls_arr)
cols = 10
rows = int(int(divmod(nbls,cols)[0] + np.ceil(divmod(nbls,cols)[1]/float(cols))))
tinfo,data,flags = C.arp.get_dict_of_uv_data(args,bls_str,opts.pol)
odata = C.sak.order_data(data,info)

if len(opts.ref) == 1: refbl = bl2tup(red_bls_arr[int(opts.ref)])
else: refbl = bl2tup(opts.ref)

fig,axarr = plt.subplots(rows,cols,sharex=True,sharey=True)
plt.suptitle('Relative to %s'%str(refbl))
for i, ax in enumerate(axarr.ravel()):
    try: bl = bl2tup(red_bls_arr[i])
    except IndexError: break 
    ax.set_title(str(bl),color='black')
    if divmod(i,cols)[-1] != 0: ax.yaxis.set_visible(False)
    if divmod(i,cols)[0] != rows-1: ax.xaxis.set_visible(False)
    if bl==refbl: continue
    
    try: ax.imshow(np.angle(odata[bl][opts.pol]*np.conj(odata[refbl][opts.pol])),vmax=np.pi,vmin=-1.*np.pi,aspect='auto',interpolation='nearest')
    except KeyError:
        print 'KeyError: %s,%s'%(bl,refbl)
        continue
fig.subplots_adjust(wspace=0,hspace=0)
plt.show()
plt.close()

"""
for i, bl in enumerate(red_bls_arr):
    ax = plt.subplots(rows,cols,i+1)
    #plt.title(str(bl)+' '+refbl,color='magenta')
    #C.sak.waterfall(odata[bl][opts.pol]*np.conj(odata[refbl][opts.pol]),\
    #                mode='phs', cmap='jet', mx=np.pi, drng=2*np.pi)
    ax.imshow(odata[bl][opts.pol]*np.conj(odata[refbl][opts.pol]),vmax=np.pi,vmin=-1.*np.pi,aspect='auto',interpolation='nearest')
    #plt.grid(0)
    if divmod(i,cols)[-1] != 0: ax.yaxis.set_visible(False)
    if divmod(i,cols)[0] != rows-1: ax.xaxis.set_visible(False)
plt.show()
"""




"""
print fengines
print N_red_bls_per_fengine
rows = len(fengines)
cols = max(N_red_bls_per_fengine)
fig = plt.figure(figsize=(16,12))
for cnt,i in enumerate(fengines):
    if opts.verb: print 'fengine %i'%i
    ants_str = C.sak.get_ants(opts.mfile,pol=opts.pol,search='f%i'%i)
    ants = map(int,ants_str.split(',')) #[0,1,2,...]
    f_bls = []
    # get all baselines on that fengine
    for j in ants:
        for k in ants:
            if j!=k:
                f_bls.append('%i_%i'%(j,k))
    # get the redundant baselines on that fengine
    pltbls = []
    for bl in f_bls:
        if bl in red_bls_arr:
            pltbls.append(bl)
    firstbl = pltbls[0]
    for pbl in pltbls:
sys.exit()
L = len(red_bls_arr)
LP = L*(L-1)/2
storage = np.zeros((L,203,19*len(sys.argv[1:])),dtype='complex128')

for b,bl in enumerate(red_bls_arr):
        t,d,f = C.arp.get_dict_of_uv_data(args,bl,opts.pol)
        od = C.sak.order_data(d,info)
        storage[b,:,:] = od
for i in range

C.sak.plot_phase_ratios(od)
plt.show()
"""


