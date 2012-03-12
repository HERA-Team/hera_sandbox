#!/usr/bin/env python
#
#  vo_cone_specfind_pca.py
#  
#
#  Created by Danny Jacobs on 3/10/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,ephem as ep
from astrogrid import ConeSearch
import vo.table as vot
from cStringIO import StringIO


o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--uri',default="ivo://CDS.VizieR/VIII/85",
    help="URI for vo data source that supports cone search. default=ivo://CDS.VizieR/VIII/85 [specfindv2]")
o.add_option('--max_degen',
    help="If there are M frequencies available, limit the covariance matrix to a size of Mx<max_degen>*M ")
opts, args = o.parse_args(sys.argv[1:])


src = a.phs.RadioFixedBody("4:13:49.93","11:11:06.1")
ep.FixedBody.compute(src,ep.J2000)
conesrch = ConeSearch(opts.uri)
votable = conesrch.execute(src.ra/ep.degree,src.dec/ep.degree,5)
votable_f = StringIO(votable)
votable_np = vot.parse_single_table(votable_f,pedantic=False)
vodata = votable_np.array

all_freqs = n.sort(n.array(list(set(vodata['nu']))))
seqs = list(set(vodata['Seq']))
print seqs
X = n.matrix(n.zeros([len(seqs),len(all_freqs)]))
for i,seq in enumerate(seqs):
    spec = n.sort(vodata[n.where(vodata['Seq']==seq)[0]],order='nu')
    flux_sc = n.log10(spec['S_nu_'])-n.average(n.log10(spec['S_nu_']))
    P = n.poly1d(n.polyfit(n.log10(spec['nu']),n.log10(spec['S_nu_']),2))
    chi_sq = n.sum(n.log10(spec['S_nu_'])-n.log10(P(n.log10(spec['nu'])))/(2*n.log10(spec['e_S_nu_'])))/(len(spec)-1)
    for j,freq in enumerate(all_freqs):
        X[i,j] = n.average(flux_sc[n.where(spec['nu']==freq)[0]]) #bin up the frequencies (accounting for degenerate measurements)
    if not opts.max_degen is None: 
        if i==len(all_freqs): break
cov = n.cov(X)
w,v = n.linalg.eig(cov)