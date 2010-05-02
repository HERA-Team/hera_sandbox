#!/usr/bin/env python
"""
A script for fitting parameters of a measurement equation given 
starting parameters in a cal file and a list of sources.  The fitter used
here is a steepest-decent filter and does not make use of priors.
"""

import aipy as a, numpy as n, sys, optparse, spead, cPickle
#import logging; logging.basicConfig(level=logging.DEBUG)

CMD_CACHE, CMD_PROCESS, CMD_FLUSH, CMD_CLEARCACHE, CMD_CLEARALL = range(5)

o = optparse.OptionParser()
o.set_usage('fitmdl_master.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True,
    cal=True, src=True, dec=True, prms=True)
o.add_option('-S', '--shared_prms', dest='shprms',
    help='Parameter listing w/ same syntax as "-P/--prms" except that all objects listed for a parameter will share an instance of that parameter.')
o.add_option('--snap', dest='snap', action='store_true',
    help='Snapshot mode.  Fits parameters separately for each integration.')
o.add_option('-q', '--quiet', dest='quiet', action='store_true',
    help='Be less verbose.')
o.add_option('--maxiter', dest='maxiter', type='float', default=-1,
    help='Maximum # of iterations to run.  Default is infinite.')
o.add_option('--xtol', dest='xtol', type='float', default=1e-10,
    help='Fractional change sought in it parameters before convergence.  Default 1e-10.')
o.add_option('--ftol', dest='ftol', type='float', default=1e-10,
    help='Fractional tolerance sought in score before convergence.  Default 1e-10.')
o.add_option('--remem', dest='remember', action='store_true',
    help='Remember values from last fit when fitting in snapshot mode.')
o.add_option('--baserx', dest='baserx', type='int', default=53000,
    help="Base port # to use for rx.  Each daemon adds it's daemon id to this to determine the actual port used for TCP transactions.")
o.add_option('--daemons', dest='daemons', default='127.0.0.1:53001',
    help='Operate in master mode, employing daemon-mode servers to do the work and collecting the results.  Should be a comma delimited list of host:port pairs to contact.')
o.add_option('--hostname', dest='hostname', default='127.0.0.1',
    help='name or IP address of the host computer running fitmdl_master')
o.add_option('--sim_autos', dest='sim_autos', action='store_true',
    help='Use auto-correlations in fitting.  Default is to use only cross-correlations.')

opts, args = o.parse_args(sys.argv[1:])

# Parse command-line options
opts.ant += ',cross'
if opts.maxiter < 0: opts.maxiter = n.Inf
# Process daemons
def parsehostport(hostport):
    host, port = hostport.split(':')
    return (host, int(port))
daemon_ip_ports = [parsehostport(w) for w in opts.daemons.split(',')]
NDAEMONS = len(daemon_ip_ports)
_igs = [spead.ItemGroup() for i in range(NDAEMONS)]
_txs = [spead.Transmitter(spead.TransportUDPtx(*ip_port)) for ip_port in daemon_ip_ports]
_igrxs = [spead.ItemGroup() for i in range(NDAEMONS)]
_iters = [spead.iterframes(spead.TransportUDPrx(opts.baserx+i)) for i in range(NDAEMONS)]

def igs(key, value, num=None):
    if num is None:
        for ig in _igs: ig[key] = value
    else: _igs[num][key] = value
def add_items(*args, **kwargs):
    num = kwargs.pop('num', None)
    if num is None:
        for ig in _igs: ig.add_item(*args, **kwargs)
    else: _igs[num].add_item(*args, **kwargs)
def send_frames(num=None):
    if num is None:
        for ig,tx in zip(_igs,_txs): tx.send_frame(ig.get_frame())
    else: _txs[num].send_frame(_igs[num].get_frame())
def ends(num=None):
    if num is None:
        for tx in _txs: tx.end()
    else: _txs[num].end()

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
NCHAN = chans.size
aa.select_chans(chans)
srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
(uvw,t,(i,j)),d = uv.read()
aa.set_jultime(t)
cat.compute(aa)

c8 = spead.mkfmt(('c',8))
u40 = spead.mkfmt(('u',40))
f64 = spead.mkfmt(('f',64))
f32 = spead.mkfmt(('f',32))
i32 = spead.mkfmt(('i',32))

add_items('command', fmt=u40, shape=[], init_val=CMD_CLEARALL)
send_frames()

add_items('baseline', fmt=spead.mkfmt(('u',20),('u',20)), shape=[])
add_items('juldate', fmt=f64, shape=[])
add_items('pol', fmt=spead.mkfmt(('u',24),('c',8),('c',8)), shape=[])
add_items('cal', fmt=c8, shape=-1, init_val=opts.cal)
add_items('src', fmt=c8, shape=-1, init_val=opts.src)
add_items('cat', fmt=c8, shape=-1, init_val=opts.cat)
add_items('respond_ip', fmt=c8, shape=-1, init_val=opts.hostname)
for i in range(NDAEMONS):
    add_items('respond_port', fmt=u40, shape=[], init_val=opts.baserx+i, num=i)
add_items('data_real', fmt=f32, shape=-1)
add_items('data_imag', fmt=f32, shape=-1)
add_items('flags', fmt=spead.mkfmt(('u',8)), shape=-1) ## Gotta be 8 b/c nflags * bits % 8 = 0
add_items('chans', fmt=i32, shape=-1) # XXX
add_items('sdf', fmt=f64, shape=[], init_val=uv['sdf']) 
add_items('sfreq', fmt=f64, shape=[], init_val=uv['sfreq']) 
add_items('nchan', fmt=u40, shape=[], init_val=uv['nchan']) 

# Prepare new data cache
# Split data across daemons by chans
chs = {}
for i in range(NDAEMONS): chs[i] = []
# Farm out chs to daemons, and give each daemon at least 1 ch
for cnt in range(max(NCHAN,NDAEMONS)):
    chs[cnt%NDAEMONS].append(chans[cnt%NCHAN])
for i in range(NDAEMONS):
    chs[i] = n.array(chs[i])
    ch = chs[i]
    print 'Assigning daemon %d/%d %d channels' % (i+1,NDAEMONS, len(ch)),
    if len(ch) == 1: print ch
    else: print
    igs('chans', ch.reshape(ch.shape + (1,)), num=i)
dmn_by_ch = {}
for d in chs:
    for c in chs[d]:
        dmn_by_ch[c] = dmn_by_ch.get(c, []) + [d]

# Figure out parameters to fit
add_items('key_list', fmt=c8, shape=-1)
add_items('prm_list', fmt=f64, shape=-1)
prms, prm_dict, shkeys = {}, {}, []
# Handle shared parameters
if opts.shprms:
    shprms = map(a.scripting.parse_prms, opts.shprms.split(','))
    for s in shprms:
        keys = s.keys(); keys.sort()
        k = keys[0]
        # Only enter one instance of shared parameters (indexed under first key)
        if prms.has_key(k): prms[k].update(s[k])
        else: prms.update({k:s[k]})
        # Add an entry to shkeys for propagating variables to other objects
        shkeys.append((keys, s[k].keys()))
# Handle normal parameters
if opts.prms:
    pd = a.scripting.parse_prms(opts.prms)
    for k in pd:
        if prms.has_key(k): prms[k].update(pd[k])
        else: prms[k] = pd[k]
for prm in prms: prm_dict[prm] = prms[prm].keys()
start_prms = aa.get_params(prm_dict)
start_prms.update(cat.get_params(prm_dict))
for obj in start_prms:
    for prm in start_prms[obj]:
        if prms[obj][prm][0] != None:
            start_prms[obj][prm] = prms[obj][prm][0]
        
prm_list, key_list = a.fit.flatten_prms(start_prms)
prm_list = n.array(prm_list)
igs('prm_list', prm_list.reshape(prm_list.shape + (1,)))
igs('key_list', cPickle.dumps(key_list))

#------------------------------------------------------------------------------------------------
# Call the optimizer

first_fit = None    # Used to normalize fit values to the starting fit
def fit_func(prm_list):
    global first_fit
    if first_fit == 0: return 0
    score, nsamples = 0, 0
    pdict = a.fit.reconstruct_prms(prm_list, key_list)
    if not opts.quiet: a.fit.print_params(pdict)
    prm_list = n.array(prm_list)
    igs('prm_list', prm_list.reshape(prm_list.shape + (1,)))
    igs('command', CMD_PROCESS)
    #print 'sending frame'
    send_frames()
    igs('command', CMD_FLUSH) # this is a hacky way to generate another frame to flush the first
    send_frames()
    for dnum,(ig_rx,iter_rx) in enumerate(zip(_igrxs,_iters)):
      for frame in iter_rx:
        #print 'got a frame back'
        ig_rx.update(frame)
        if ig_rx['valid']: break
      #print 'Got response from daemon %d/%d' % (dnum+1, NDAEMONS), ig_rx['score'], ig_rx['nsamples']
      score += ig_rx['score']
      nsamples += ig_rx['nsamples']
    #print 'Total:', score, nsamples
    #if not opts.quiet: print score, nsamples
    if nsamples == 0:
        first_fit = 0.
        return 0.
    score = n.sqrt(score / nsamples)
    if first_fit is None: first_fit = score
    if not opts.quiet:
        print
        print 'Score:', score, 
        print '(%2.2f%% of %f)' % (100 * score / first_fit, first_fit)
        print '-' * 70
    return score / first_fit

if not opts.snap:
    times = []
    for uvfile in args:
        print 'Sending data from', uvfile
        uv = a.miriad.UV(uvfile)
        a.scripting.uv_selector(uv, opts.ant, opts.pol)
        uv.select('decimate', opts.decimate, opts.decphs)
        for (uvw,t,(i,j)),d,f in uv.all(raw=True):
            if len(times) == 0 or times[-1] != t: times.append(t)
            active_daemons = {}
            for dnum in range(NDAEMONS):
                ch = chs[dnum]
                # If only 1 ch per daemon, then we start dividing data by time
                if len(ch) == 1:
                    dmns = dmn_by_ch[ch[0]]
                    # If it's not this guy's turn for data, get out
                    if dnum != dmns[len(times) % len(dmns)]: continue
                active_daemons[dnum] = None
                _d = d.take(ch)
                _d = _d.reshape(_d.shape + (1,))
                _f = f.take(ch)
                #print 'Sending', dnum, ch, t, (i,j)
                igs('data_real', _d.real, num=dnum)
                igs('data_imag', _d.imag, num=dnum)
                igs('flags', _f.reshape(_f.shape + (1,)), num=dnum)
            igs('pol', (0,) + tuple(a.miriad.pol2str[uv['pol']]))
            igs('juldate', t)
            igs('baseline', (i,j))
            igs('command', CMD_CACHE)
            for dnum in active_daemons.keys(): send_frames(dnum)
    rv = a.optimize.fmin(
        fit_func, prm_list,
        full_output=1, disp=0,
        maxfun=opts.maxiter, maxiter=n.Inf, 
        ftol=opts.ftol, xtol=opts.xtol
    )
    prms,score = rv[:2]
    prms = a.fit.reconstruct_prms(prms, key_list)
    print
    a.fit.print_params(prms)
    print 'Score:', score * first_fit, 
    print '(%2.2f%% of %f)' % (100 * score, first_fit)
    print '------------------------------------------------------------'
else:
    curtime, cnt = None, 0
    for uvfile in args:
        uv = a.miriad.UV(uvfile)
        a.scripting.uv_selector(uv, opts.ant, opts.pol)
        uv.select('decimate', opts.decimate, opts.decphs)
        for (uvw,t,(i,j)),d,f in uv.all(raw=True):
            if curtime != t:
                if curtime != None:
                    print 'Time:', curtime
                    print 'Iter:', cnt
                    first_fit = None
                    rv = a.optimize.fmin(
                        fit_func, prm_list,
                        full_output=1, disp=0,
                        maxfun=opts.maxiter, maxiter=n.Inf, 
                        ftol=opts.ftol, xtol=opts.xtol
                    )
                    prms,score = rv[:2]
                    prms = a.fit.reconstruct_prms(prms, key_list)
                    print
                    a.fit.print_params(prms)
                    print 'Score:', score * first_fit, 
                    print '(%2.2f%% of %f)' % (100 * score, first_fit)
                    print '-' * 70
                    if opts.remember:
                        prm_list, key_list = a.fit.flatten_prms(prms)
                    igs('command', CMD_CLEARCACHE)
                    send_frames()
                curtime, cnt = t, cnt+1
            # Cache the data for the next time step
            for dnum in range(NDAEMONS):
                ch = chs[dnum]
                _d = d.take(ch)
                _d = _d.reshape(_d.shape + (1,))
                _f = f.take(ch)
                igs('data_real', _d.real, num=dnum)
                igs('data_imag', _d.imag, num=dnum)
                igs('flags', _f.reshape(_f.shape + (1,)), num=dnum)
            igs('pol', (0,) + tuple(a.miriad.pol2str[uv['pol']]))
            igs('juldate', t)
            igs('baseline', (i,j))
            igs('command', CMD_CACHE)
            send_frames()

ends()
