#! /usr/bin/env python

import aipy as a, numpy as n, spead, optparse, sys, cPickle
#import logging; logging.basicConfig(level=logging.DEBUG)

CMD_CACHE, CMD_PROCESS, CMD_FLUSH, CMD_CLEARCACHE, CMD_CLEARALL = range(5)

o = optparse.OptionParser()
o.add_option('-p','--port', type='int', help='UDP port to listen to.')
opts, args = o.parse_args(sys.argv[1:])

ig = spead.ItemGroup()
dbuf = {}
aa, cat, tx, ig_tx = None, None, None, None

f64 = spead.mkfmt(('f',64))
i40 = spead.mkfmt(('i',40))

while True:
    tport = spead.TransportUDPrx(opts.port, pkt_count=1024)
    print 'waiting for a connection'
    for frame in spead.iterframes(tport):
        ig.update(frame)
        if ig['command'] == CMD_FLUSH:
            #print 'flushing'
            continue
        if ig['command'] == CMD_CLEARALL:
            print 'clearing all'
            dbuf, aa, cat, ig_tx = {}, None, None, None
        if ig['command'] == CMD_CLEARCACHE:
            print 'clearing cache'
            dbuf = {}
        elif ig['command'] == CMD_CACHE:
            print 'caching data'
            cal = ig['cal'].tostring()
            if aa is None:
                exec('import %s; reload(%s)' % (cal, cal)) # hack to force reload of a changed cal file
                aa = a.cal.get_aa(cal, ig['sdf'], ig['sfreq'], ig['nchan'])
                print ig['chans']
                aa.select_chans(ig['chans'][:,0])
                aa.set_jultime(ig['juldate'])
            if cat is None:
                exec('import %s; reload(%s)' % (cal, cal)) # hack to force reload of a changed cal file
                srclist, cutoff, catalogs = a.scripting.parse_srcs(ig['src'].tostring(), ig['cat'].tostring())
                cat = a.cal.get_catalog(cal, srclist, cutoff, catalogs)
                cat.compute(aa)
                mfq = cat.get('mfreq')
            if ig_tx is None:
                tport_tx = spead.TransportUDPtx(ig['respond_ip'].tostring(), ig['respond_port'])
                tx = spead.Transmitter(tport_tx)
                ig_tx = spead.ItemGroup()
                ig_tx.add_item('score', fmt=f64)
                ig_tx.add_item('nsamples', fmt=f64)
                ig_tx.add_item('valid', fmt=i40, init_val=0)
            # Finally, cache this data
            t = ig['juldate']
            bl = a.miriad.ij2bl(*ig['baseline'])
            pol = ''.join(ig['pol'][-2:])
            d = ig['data_real'].squeeze() + 1j * ig['data_imag'].squeeze()
            f = ig['flags'].squeeze()
            print t, a.miriad.bl2ij(bl)#, d
            if not dbuf.has_key(t): dbuf[t] = {}
            if not dbuf[t].has_key(pol): dbuf[t][pol] = {}
            dbuf[t][pol][bl] = (d, f, n.where(f, 0, n.abs(d)**2).sum())
        elif ig['command'] == CMD_PROCESS:
            #print 'unpacking parameters'
            # Get parameters for this iteration
            key_list = cPickle.loads(ig['key_list'].tostring())
            prms = a.fit.reconstruct_prms(ig['prm_list'], key_list)
            #print prms
            aa.set_params(prms)
            cat.set_params(prms)
            a1,a2,th = cat.get('srcshape')
            # Process data from cache
            #print 'simulating model'
            score,nsamples = 0, 0
            for t in dbuf:
                aa.set_jultime(t)
                cat.compute(aa)
                eqs = cat.get_crds('eq', ncrd=3)
                flx = cat.get_jys()
                dra,ddec = cat.get('ionref')
                aa.sim_cache(eqs, flx, mfreqs=mfq, ionrefs=(dra,ddec), srcshapes=(a1,a2,th))
                for pol in dbuf[t]:
                  for bl in dbuf[t][pol]:
                    i,j = a.miriad.bl2ij(bl)
                    d,f,nsamp = dbuf[t][pol][bl]
                    sim_d = aa.sim(i, j, pol=pol)
                    difsq = n.abs(d - sim_d)**2
                    difsq = n.where(f, 0, difsq)
                    score += difsq.sum()
                    nsamples += nsamp
            print 'returning result to (%s:%d)' % (ig['respond_ip'].tostring(), ig['respond_port'])
            # Return the result
            ig_tx['score'] = score
            ig_tx['nsamples'] = nsamples
            ig_tx['valid'] = 1
            tx.send_frame(ig_tx.get_frame())
            ig_tx['valid'] = 0
            tx.send_frame(ig_tx.get_frame())
