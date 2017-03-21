#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse, re, os, random

NBOOT = 400
MEDIAN = True
CLIP = False
LO,HI = 40,320
#LO,HI = 40,600
args = sys.argv[1:]

pk_vs_t = {}
err_vs_t = {}
temp_noise_var = {}
nocov_vs_t = {}
for filename in args:
    print 'Reading', filename
    f = n.load(filename)
    kpl,cmd = f['kpl'], f['cmd']
    path = os.path.dirname(filename)
    if not pk_vs_t.has_key(path):
        print '   ', path
        print '   ', cmd
        pk_vs_t[path] = []
        err_vs_t[path] = []
        temp_noise_var[path] = []
        nocov_vs_t[path] = []
    pk_vs_t[path].append(f['pk_vs_t'])
    scalar = f['scalar']
    nocov_vs_t[path].append(f['nocov_vs_t'])

paths = pk_vs_t.keys()
k0 = n.abs(kpl).argmin()

pk_2d = n.array([pk_vs_t[path] for path in paths]) # (bltype,bootstraps,kpls,times)
#pk_2d = n.concatenate([pk_2d[...,11:,:], pk_2d[...,9::-1,:]], axis=-1)
pk_2d = n.transpose(pk_2d, [0,3,1,2]) # (bltype,times,boot,kpls)
#pk_2d = pk_2d[:,::8,...].copy()
pk_2d = pk_2d[:,::2,...].copy()
pk_2d.shape = (pk_2d.shape[0]*pk_2d.shape[1],) + pk_2d.shape[2:] # (modes, boot, kpls)
print pk_2d.shape

acc = n.concatenate([[1], n.around(2**n.arange(1,n.log2(pk_2d.shape[0]), .5))])
avg_var = n.zeros((acc.size,)+pk_2d.shape[1:], dtype=n.float)
med_var = n.zeros((acc.size,)+pk_2d.shape[1:], dtype=n.float)
print avg_var.shape
SHUFFLES = 200
print acc
for cnt in xrange(SHUFFLES):
    print cnt
    for i,t in enumerate(acc):
        avg_var[i] += n.abs(n.average(pk_2d[:t], axis=0).real)**2
        med_var[i] += n.abs( n.median(pk_2d[:t], axis=0).real)**2
    n.random.shuffle(pk_2d)
#avg_var = n.average(avg_var, axis=1) / SHUFFLES # avg over boots
#med_var = n.average(med_var, axis=1) / SHUFFLES # avg over boots
avg_var /= SHUFFLES
med_var /= SHUFFLES


#cmap = p.get_cmap('jet')
cmap = p.get_cmap('brg')
#for cnt, kpl in enumerate(xrange(1,pk_2d.shape[-1])):
for cnt in xrange(0,pk_2d.shape[-1]):
    k = n.abs(kpl[cnt])
    #avg_var[:,kpl] /= avg_var[0,kpl]
    #med_var[:,kpl] /= med_var[0,kpl]
    if cnt in [9, 10, 11]: continue
    print k
    #c = colors[cnt % len(colors)]
    c = cmap((n.abs(cnt-10)-2)/8.)
    #p.loglog(acc,k**3/(2*n.pi**2)*n.sqrt(avg_var[:,cnt]), '--', color=c)
    #p.loglog(acc,k**3/(2*n.pi**2)*n.sqrt(med_var[:,cnt]), '-', color=c)
    for b in xrange(pk_2d.shape[1]):
        p.subplot(121); p.loglog(acc,k**3/(2*n.pi**2)*n.sqrt(avg_var[:,b,cnt]), '-', alpha=.2, color=c)
        p.subplot(122); p.loglog(acc,k**3/(2*n.pi**2)*n.sqrt(med_var[:,b,cnt]), '-', alpha=.2, color=c)

p.subplot(121); p.xlim(1,1e3); p.ylim(1e1,1e5); p.grid()
p.title('Mean')
p.ylabel(r'$\langle k^3/2\pi^2\ P(k)\rangle\ [{\rm mK}^2]$', fontsize=14)
p.xlabel('Modes integrated')

p.subplot(122); p.xlim(1,1e3); p.ylim(1e1,1e5); p.grid()
p.title('Median')
p.xlabel('Modes integrated')
p.setp(p.gca().get_yticklabels(), visible=False)
p.show()
