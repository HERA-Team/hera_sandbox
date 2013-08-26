#! /usr/bin/env python
"""
freq_average_beamform.py

Takes a set of npz files and produces Nchan sets of averaged files
"""
import aipy as a, numpy as n
import optparse, sys, re



o = optparse.OptionParser()
o.set_description(__doc__)
o.add_option('--chan_width',default=10.,
    help="channel width in MHz")
o.add_option('--err_wgt',action='store_true',
    help='weight by the inverse of the coarse channel variance')
opts,args = o.parse_args(sys.argv[1:])    


for file in args:
    print file,
    D = n.load(file)
    freq = D['freq']*1e3

    #use the nearest integeger number of channels
    Ncoarse = int(n.ceil((freq.max() - freq.min())/opts.chan_width))
    Nfine = freq.size/Ncoarse
    df = (freq.max() - freq.min())/Ncoarse
    Nedge = freq.size - Ncoarse*Nfine
    #cut off the "remainder" channels after choosing an integer number of bins
    chans = n.arange(Nedge/2,Ncoarse*Nfine+1)
    #grab the spectrum and average it down
    SPECTRUM = n.ma.masked_where(D['flags'],D['spec'])[:,chans]
    SPECTRUM.shape = (SPECTRUM.shape[0],Ncoarse,Nfine)
    AVG_SPECTRUM = n.ma.mean(SPECTRUM,axis=1)
    ERR_SPECTRUM = n.ma.std(SPECTRUM,axis=1)

    #output the results
    MFREQs = n.arange(Ncoarse)*df + freq[chans].min() + df/2.
    outD = {}
    for thing in D:
        if not thing in ['spec','wgts']:
            outD[thing] = D[thing]

    for i,f in enumerate(MFREQs):
        outfile = "%s_%5.1fMHz"%(file[:-4],f)
        print outfile,
        if opts.err_wgt:
            weights = D['wgts']/ERR_SPECTRUM[:,i]
        else:
            weights = D['wgts']
        n.savez(outfile,spec=AVG_SPECTRUM[:,i].filled(0),wgt=weights.filled(0),**outD)
    print 

