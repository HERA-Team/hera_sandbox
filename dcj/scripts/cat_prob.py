#!/usr/bin/env python
#
#  cat_prob.py
#  
#
#  Created by Danny Jacobs on 4/23/10.
#  PAPER Project
#

import aipy as a, pylab as p,math as m
import sys, optparse,vo.table as vot,atpy
from numpy import *
from pylab import *
import numpy as n
from astrogrid import ConeSearch
from cStringIO import StringIO


o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])


dNds = n.loadtxt('nvss_dnds.txt')
mask = n.zeros_like(dNds)
mask[n.argwhere(dNds==0)] = 1
dNds = n.ma.array(dNds,mask=mask)
fit = n.ma.polyfit(dNds[:15,0],n.log10(dNds[:15,1]),1)
E = lambda s: 10**fit[1]*s**fit[0]
G = lambda s,ds: (10**fit[1]/fit[0])*((s+ds)**(fit[0]+1) - s**(fit[0]+1))
figure(34)
clf()
loglog(10**dNds[:,0],dNds[:,1],'.')
loglog(10**dNds[:15,0],dNds[:15,1],'+')
loglog(10**dNds[:,0],E(10**dNds[:,0]),'k')

#fit = array([-2.64027037,  4.71238127])
deg = pi/180
#dNdO = lambda s: 10**fit[1]*s**(fit[0]-1)*0.1
As = 7.2
s = logspace(-1,2)
s0 = 5
ds0 = 1.2
r0 = 10.0*deg
r = logspace(-1,2)
dr = 1.0*deg
cbeam = (0.01*deg)
pbeam = (0.6*deg)
fds = 0.1
#possible candidates
#10068
#srcs  = array([(28.6,2.7,0.078*a.const.arcsec,4.02*a.const.arcsec),
#            (0.75,0.04,0.1*a.const.arcsec,0.1*a.const.arcsec),
#            (3.1,0.8,0.0,1e-5*deg),
#            (4.2,2.2,0.0,1e-5*deg),
#            (8.1,4.0,0.0,1e-5*deg)],dtype({'names':['S','e_S','r','dr'],
#            'formats':(float,)*4}))
#srcs['r'] = n.random.uniform(0,0.5*deg,size=len(srcs))

R,S = meshgrid(r,s)
def Prob(L): return (1-exp(-L))
#def Prob2(L): return L**2*exp(-L)/2
print "number of sources within %f%% of %f Jy within %f deg:"%(fds*100,s0,r0),
print (pi*r0**2)*G(s0,fds*s0)
print "given a probability of finding that source",
print Prob((pi*r0**2)*G(s0,fds*s0))
#print "probability of finding two sources within the same flux and position range",
#print Prob2((r0*dr)*G(s0,fds*s0))
#print 1-E(s0)*As*2*pi*dr*0.1*pi/180.*s0*fds*exp(-E(s0)*As*2*pi*dr*0.1*pi/180.*s0*fds)
#Ps = (1 - As*2*pi*dr*R*pi/180.*S*fds*E(S)*exp(-As*2*pi*dr*R*pi/180.*S*fds*E(S)))
for s0 in [1,3.5,5,10,19,35,50,100]:
    
    P = (1-Prob(pi*(R*deg)**2*G(S,S*fds)))*(1-Prob(pi*(R*deg)**2.*G(s0,s0*fds)))
    P2   = exp(-(log10(s0)-log10(S))**2/(2*((ds0/s0)**2+(fds)**2)))
    P3    =exp(-(R*pi/180)**2/(2*(dr**2+cbeam**2)))
    fig = figure(56)
    clf()
    conf = 0.75
    cs1 = contour(log10(R),log10(S),P,[conf])#1-array([0.25,0.5,0.75,0.95,0.99]))
    cs2 = contour(log10(R),log10(S),P2,[conf])
    cs3 = contour(log10(R),log10(S),P3,[conf])
    contour(log10(R),log10(S),P2*P3,[conf])
    cs = contour(log10(R),log10(S),P*P2*P3,[conf],linewidths=3)
    xlabel('log radius [deg]') 
    ylabel('log flux [jy]')
    draw()
    fig.savefig('cat_prob_P_'+str(s0)+'_'+str(s0*fds)+'.png')

votable_np = vot.parse_single_table(open('../../cals/test.vot'),pedantic=False)
VO = atpy.Table('../../cals/test.vot')
VO.add_empty_column('P_dnds',dtype=float,unit='prob',description='prior prob given number counts')
VO.add_empty_column('P_r',dtype=float,unit='prob',description='radial id probability')
VO.add_empty_column('P_S',dtype=float,unit='prob',description='flux id probability')
VO.add_empty_column('P',dtype=float,unit='prob',description='total probability = P_dnds*P_r*P_S')



for exseq in list(set(votable_np.array['PAPER_seq'])):
    data = votable_np.array[n.where(votable_np.array['PAPER_seq']==exseq)[0]]
    s0 =n.max(data['S_nu_paper'])
    ds0 = n.max(data['e_S_nu_paper'])
    seqs = list(set(data['Seq']))
    print "looking for source of flux %7.3f +/- %7.3f"%(s0,ds0)
    print '\t'.join(('seq','S','dS','R','dR','P_dnds','P_S','P_r','Ptot'))
    for seq in seqs:
        S = n.max(data[n.where(data['Seq']==seq)[0]]['S_sf'])#max is just a trick to get one element of an array where all the elements should be the same
        dS = n.max(data[n.where(data['Seq']==seq)[0]]['e_S_sf'])
        R = n.average(data[n.where(data['Seq']==seq)[0]]['_r'])
        dR = n.average(data[n.where(data['Seq']==seq)[0]]['beam'])*a.const.arcsec
        P_dnds = (1-Prob(pi*(R*pi/180.)**2*G(S,S*fds)))*(1-Prob(pi*(R*pi/180)**2.*G(s0,ds0)))
        P_strict = 1
        P_r   = exp(-(log10(s0)-log10(S))**2/(2*((ds0/s0)**2+(dS/S)**2)))/sqrt(2*pi*((ds0/s0)**2+(dS/S)**2))*sqrt((ds0/s0)**2+(dS/S)**2)
        P_S    =exp(-(R*deg)**2/(2*(dr**2+dR**2)))/sqrt(2*pi*(dr**2+dR**2))*sqrt(dr**2+dR**2)
        print '\t'.join(('%d','%4.2f','%4.2f','%6.4f','%4.3f','%4.3f','%3.2e','%3.2e','%3.2e'))%(seq,S,dS,R,dR,P_dnds,P_r,P_S,P_dnds*P_r*P_S)    
#        VO.data['P_dnds'][n.where(VO.data['Seq']==seq)[0]] = P_dnds
#        VO.data['P_r'][n.where(VO.data['Seq']==seq)[0]] = P_r
#        VO.data['P_S'][n.where(VO.data['Seq']==seq)[0]] = P_S
#        VO.data['P'][n.where(VO.data['Seq']==seq)[0]] = P_dnds*P_r*P_S
#        
#VO.write('../../cals/test_cat_prob.vot')
