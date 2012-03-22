#!/usr/bin/env python
"""
A script for fitting parameters of a measurement equation given 
starting parameters in a cal file and a list of sources.  The fitter used
here is a steepest-decent filter and does not make use of priors.
"""

import aipy as a, numpy as n, sys, os, optparse

o = optparse.OptionParser()
o.set_usage('fitimg.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True,
    cal=True, dec=True, prms=True)
o.add_option('-i', '--imfits', dest='imfits',
    help='FITS image file to use to model visibilities')
o.add_option('-b', '--bmfits', dest='bmfits', default=None,
    help='FITS beam file to use to model visibilities')
o.add_option('-S', '--shared_prms', dest='shprms',
    help='Parameter listing w/ same syntax as "-P/--prms" except that all objects listed for a parameter will share an instance of that parameter.')
o.add_option('--snap', dest='snap', action='store_true',
    help='Snapshot mode.  Fits parameters separately for each integration.')
o.add_option('-q', '--quiet', dest='quiet', action='store_true',
    help='Be less verbose.')
o.add_option('--no_w', dest='no_w', action='store_true',
    help="Don't use W projection.")
o.add_option('--mfreq', dest='freq', type='float', default=.150,
    help="Frequency at which FITS image was made.")
o.add_option('--maxiter', dest='maxiter', type='float', default=-1,
    help='Maximum # of iterations to run.  Default is infinite.')
o.add_option('--xtol', dest='xtol', type='float', default=1e-10,
    help='Fractional change sought in it parameters before convergence.  Default 1e-10.')
o.add_option('--ftol', dest='ftol', type='float', default=1e-10,
    help='Fractional tolerance sought in score before convergence.  Default 1e-10.')
o.add_option('--remem', dest='remember', action='store_true',
    help='Remember values from last fit when fitting in snapshot mode.')
o.add_option('--baseport', dest='baseport', type='int', default=53000,
    help="Base port # to use for tx/rx.  Each daemon adds it's daemon id to this to determine the actual port used for TCP transactions.")
o.add_option('--daemon', dest='daemon', type='int', 
    help='Operate in daemon mode, opening a TCP Server to handle requests on the specified increment to the base port.')
o.add_option('--master', dest='master', 
    help='Operate in master mode, employing daemon-mode servers to do the work and collecting the results.  Should be a comma delimited list of host:daemonid pairs to contact.  Daemon ID will be added to baseport to determine actual port used for TCP transactions.')
o.add_option('--sim_autos', dest='sim_autos', action='store_true',
    help='Use auto-correlations in fitting.  Default is to use only cross-correlations.')

opts, args = o.parse_args(sys.argv[1:])
dim, kwds = a.img.from_fits(opts.imfits); dim = dim.squeeze()
if not opts.bmfits is None:
    dbm, kwds = a.img.from_fits(opts.bmfits); dbm = dbm.squeeze()
else:
    dbm = n.zeros_like(dim)
    dbm[dbm.shape[0]/2,dbm.shape[1]/2] = 1
if True:
    DIM = dim.shape[0]
    DIM2 = DIM/2
    new_dim = n.zeros((DIM*2,DIM*2), dtype=dim.dtype)
    new_dim[DIM-DIM2:DIM+DIM2,DIM-DIM2:DIM+DIM2] = dim
    dim = a.img.recenter(new_dim, (DIM,DIM))
    new_dbm = n.zeros((DIM*2,DIM*2), dtype=dbm.dtype)
    new_dbm[DIM-DIM2:DIM+DIM2,DIM-DIM2:DIM+DIM2] = dbm
    dbm = a.img.recenter(new_dbm, (DIM,DIM))
else:
    dim = a.img.recenter(dim, (dim.shape[0]/2,dim.shape[0]/2))
    dbm = a.img.recenter(dbm, (dbm.shape[0]/2,dbm.shape[0]/2))

src = a.phs.RadioFixedBody(kwds['ra']*a.img.deg2rad, kwds['dec']*a.img.deg2rad)
if opts.no_w: Img = a.img.Img
else: Img = a.img.ImgW
size = 1 / (kwds['d_ra'] * a.img.deg2rad)
res = size / dim.shape[0]
print 'SIZE=%f, RES=%f' % (size, res)
im = Img(size, res, mf_order=0)
assert(im.uv.shape == dim.shape)
im.uv = n.fft.fft2(dim).astype(n.complex64)
im.bm = [n.fft.fft2(dbm).astype(n.complex64)]

# Parse command-line options
uv = a.miriad.UV(args[0])
opts.ant += ',cross'
a.scripting.uv_selector(uv, opts.ant, opts.pol)
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
aa.select_chans(chans)
(uvw,t,(i,j)),d = uv.read()
aa.set_jultime(t)
src.compute(aa)
del(uv)
if opts.maxiter < 0: opts.maxiter = n.Inf

# Figure out parameters to fit
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
for obj in start_prms:
    for prm in start_prms[obj]:
        if prms[obj][prm][0] != None:
            start_prms[obj][prm] = prms[obj][prm][0]
        
prm_list, key_list = a.fit.flatten_prms(start_prms)

first_fit = None    # Used to normalize fit values to the starting fit
dbuf = None

# The function to be optimized
def fit_func(prms, filelist, decimate, decphs):
    global first_fit, dbuf
    if first_fit == 0: return 0
    prms = a.fit.reconstruct_prms(prms, key_list)
    # Propagate shared params
    for (skey,sprm) in shkeys:
        k = skey[0]
        for k2 in skey[1:]:
            if not prms.has_key(k2): prms[k2] = {}
            for sp in sprm:
                prms[k2][sp] = prms[k][sp]
    if not opts.quiet: a.fit.print_params(prms)
    aa.set_params(prms)
    score,nsamples = 0.,0.
    # Cache data from file to avoid hitting disk excessively
    if dbuf is None:
        if not opts.quiet: print 'Caching data...'
        dbuf = {}
        for uvfile in filelist:
            sys.stdout.write('.') ; sys.stdout.flush()
            uv = a.miriad.UV(uvfile)
            a.scripting.uv_selector(uv, opts.ant, opts.pol)
            uv.select('decimate', decimate, decphs)
            for (uvw,t,(i,j)),d,f in uv.all(raw=True):
                if not dbuf.has_key(t): dbuf[t] = {}
                if not opts.sim_autos and i == j: continue
                bl = a.miriad.ij2bl(i,j)
                d = d.take(chans)
                f = f.take(chans)
                dbuf[t][bl] = (d, f, 
                        n.where(f, 0, n.abs(d)**2).sum(),
                        a.miriad.pol2str[uv['pol']])
        if not opts.quiet:
            samp, vsamp, wgts = 0, 0, 0.
            for t in dbuf:
              for bl in dbuf[t]:
                samp += len(dbuf[t][bl][1])
                vsamp += n.logical_not(dbuf[t][bl][1]).astype(n.int).sum()
                wgts += dbuf[t][bl][2]
            print 'Cache summary:'
            print '   %d samples' % samp
            print '   %d valid' % vsamp
            print '   %f sum weights' %  wgts
            sys.stdout.flush()
    # Process data from cache
    for t in dbuf:
        aa.set_jultime(t)
        src.compute(aa)
        for bl in dbuf[t]:
            i,j = a.miriad.bl2ij(bl)
            u,v,w = aa.gen_uvw(i, j, src=src)
            d,f,nsamp,pol = dbuf[t][bl]
            if nsamp == 0: continue
            sim_d = im.get((u,v,w)) * n.conj(aa.gen_phs(src, i, j))
            if True: passband = aa.passband(i,j)
            else:
                afreqs = aa.get_afreqs()
                valid = n.where(n.abs(d) > 0, 1, 0)
                passband = n.polyfit(n.log10(afreqs.compress(valid)), n.log10(n.abs(d.compress(valid)) / n.abs(sim_d.compress(valid))), deg=3)
                print passband
                passband = 10**n.polyval(passband, n.log10(afreqs))
            sim_d *= passband
            if False:
                import pylab as P
                P.subplot(131)
                P.imshow(n.log10(n.abs(im.uv) / n.abs(im.bm)).squeeze(), vmax=4, vmin=2)
                P.subplot(132)
                P.semilogy(n.abs(d),'k')
                P.semilogy(n.abs(sim_d),'r')
                P.subplot(133)
                P.plot(n.angle(d), 'k')
                P.plot(n.angle(sim_d), 'r')
                P.plot(n.angle(d * n.conj(sim_d) / n.abs(sim_d)), 'g')
                #P.plot(n.real(d - sim_d),'g')
                P.show()
            difsq = n.abs(d - sim_d)**2
            difsq = n.where(f, 0, difsq)
            score += difsq.sum()
            nsamples += nsamp
    if opts.daemon: return score, nsamples
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

# Call the optimizer
if opts.daemon:
    import SocketServer, struct
    class TCPServer(SocketServer.TCPServer):
        allow_reuse_address = True
    class FitHandler(SocketServer.BaseRequestHandler):
        def setup(self): print self.client_address, 'connected'
        def finish(self): print self.client_address, 'disconnected'
        def handle(self):
            data = self.request.recv(struct.calcsize('d')*len(prm_list))
            data = struct.unpack('<%dd' % (len(prm_list)), data)
            score, nsamples = fit_func(data, args, opts.decimate, opts.decphs)
            print 'Returning score=%f, nsamples=%f' % (score, nsamples)
            rv = struct.pack('<dd', score, nsamples)
            self.request.send(rv)
            sys.stdout.flush()
    s = TCPServer(('', opts.baseport + opts.daemon), FitHandler)
    print 'Starting daemon on TCP port %d' % (opts.baseport + opts.daemon)
    sys.stdout.flush()
    s.serve_forever()
elif opts.master:
    import socket, struct
    def parsehostport(hostport):
        host, port = hostport.split(':')
        return (host, opts.baseport + int(port))
    hostports = [parsehostport(w) for w in opts.master.split(',')]
    def fit_func(prms):
        global first_fit
        if first_fit == 0: return 0
        pdict = a.fit.reconstruct_prms(prms, key_list)
        if not opts.quiet: a.fit.print_params(pdict)
        data = struct.pack('<%dd' % (len(prms)), *prms)
        socks = []
        for hostport in hostports:
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.connect(hostport)
            sock.send(data)
            socks.append(sock)
        score, nsamples = 0, 0
        for sock in socks:
            data = sock.recv(1024)
            scr, nsp = struct.unpack('<dd', data)
            score += scr; nsamples += nsp
            if not opts.quiet: print score, nsamples
        score = n.sqrt(score / nsamples)
        if nsamples == 0:
            first_fit = 0.
            return 0.
        if first_fit is None: first_fit = score
        if not opts.quiet:
            print
            print 'Score:', score, 
            print '(%2.2f%% of %f)' % (100 * score / first_fit, first_fit)
            print '-' * 70
        return score / first_fit
    print 'Starting in master mode...'
    rv = a.optimize.fmin(
        fit_func, prm_list,
        #args=(args, opts.decimate, opts.decphs),
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
            
    

elif not opts.snap:
    rv = a.optimize.fmin(
        fit_func, prm_list,
        args=(args, opts.decimate, opts.decphs),
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
    for uvfile in args:
        # Figure out what times are in the file
        uv = a.miriad.UV(uvfile)
        a.scripting.uv_selector(uv, opts.ant, opts.pol)
        uv.select('decimate', opts.decimate, opts.decphs)
        times = [0]
        for (uvw,t,(i,j)),d,f in uv.all(raw=True):
            if times[-1] != t: times.append(t)
        times = times[1:]
        del(uv)
        decimate = len(times) * opts.decimate
        # Fit each time separately
        for cnt, t in enumerate(times):
            print 'Time:', t
            print 'Iter: %d / %d' % (cnt+1, len(times))
            first_fit,dbuf = None,None
            rv = a.optimize.fmin(
                fit_func, prm_list,
                args=([uvfile], decimate, opts.decimate*cnt + opts.decphs),
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
            if opts.remember: prm_list, key_list = a.fit.flatten_prms(prms)
