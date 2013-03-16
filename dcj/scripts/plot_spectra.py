#!/usr/bin/env python
#
#  img_cat_flux.py
#  
#
#  Created by Danny Jacobs on 10 April 2012.
#  PAPER Project
#


import aipy as a, numpy as n, math as m, ephem, sys,optparse, pyfits as pf
from pylab import *
import atpy


o = optparse.OptionParser()
#a.scripting.add_standard_options(o, cal=True,src=True)
#o.add_option('-r',dest='radius',default=1.,type='float',
#    help="Analysis region around each input source, in degrees. [1]")
#o.add_option('-o',dest='outfile',
#    help="Output the catalog in a 'standard format' [None]")

opts, files = o.parse_args(sys.argv[1:])


n.set_printoptions(precision=2,suppress=True)

#READ in the data
mytables = {}
for F in files:
    mytables[F] = atpy.Table()
    mytables[F].read(F)

#for each source get a table of data

#a set of all the source names in the files
srcnames = list(set(n.hstack([mytables[F]['Name'] for F in files])))

#go through each table and grab the data for each source
mydata = {}
for src in srcnames:
    for F in files:
        t = mytables[F] #get the current table
        try:
            mydata[src].append(t.where(t['Name']==src))
        except(KeyError):
            mydata[src] = t.where(t['Name']==src)

#make a nice subplot
nsrc = len(srcnames)
R = min(6,n.ceil(n.sqrt(nsrc)))
C = min(6,n.ceil(nsrc/R))
figcount=0
for i,src in enumerate(sorted(srcnames)):
    j = i%(R*C)+1
    if j==1: 
        figure(figsize=(20,10))
        print "making new figure"
        sys.stdout.flush()
        figcount+=1
    ax = subplot(R,C,j)
    freqs = mydata[src]['freq']/1e6
    semilogy(freqs,mydata[src]['S_nu_'],'.k')
    semilogy(freqs,mydata[src]['S_cat'],'--k')
    ylim([1,100])
    xlabel('Freq [MHz]')
    ylabel('Flux [JY]')
    text(0.1,0.8,src,transform=ax.transAxes)
    if ( (j-1) % C):yticks([]);ylabel('')
    if (R*C-C)>(j): xticks([]);xlabel('')
    if j==R*C:
        subplots_adjust(hspace=0,wspace=0)
        draw()
        savefig('plot_spectra_n%d.png'%figcount)
        print "saving new figure"
        sys.stdout.flush()


