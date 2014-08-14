#!/bin/sh
''''exec python -u -- "$0" ${1+"$@"} # '''
# vi: syntax=python
import numpy as n, pylab as p
import aipy as a
import sys
import optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, src=True)
o.add_option('--lst', dest='lst', action='store_true',
    help='plot against lst, in addtion to jd')
opts, args = o.parse_args(sys.argv[1:])

if not opts.src is None:
    srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
else: srclist,cutoff,catalogs = [], None, []; src=None

cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
aa = a.cal.get_aa(opts.cal, n.array([.15]))
if not opts.src is None:
    src = cat[opts.src]
#    print 'Avoiding src : ', src


#PLOT_GAIN = False
PLOT_GAIN = 1
colors = 'kkbgcmry'

dly = {}
gain = {}
times = []
lsts = []
days = []
jday =0 
for filename in args[1:]:
    print 'Reading', filename
    f = open(filename)
    curday = filename.split('.')[1]
    if curday != jday:
        dly[curday] = {}
        gain[curday] = {}
        days.append(curday)
    npz = n.load(f)
    C_phs = npz['C_phs']
    C_amp = 10**npz['C_amp']
    antpos = npz['antpos']
    time = npz['time']
#    print time[0], time[-1]
#    pols = npz['pols']
    pols = ['xx']
    aa.set_jultime(time)
    cat.compute(aa)
    if not src is None:
    #    print src.alt
        if src.alt>0 : 
            print 'Alt', src.alt
            continue
    lsts.append(aa.sidereal_time()*(12./n.pi))
    times.append(time)
    try:
        print times[-1] - times[-2]
    except: pass
    for pi, pol in enumerate(pols):
        if not dly[curday].has_key(pi): dly[curday][pi],gain[curday][pi] = {},{}
        for i,tau,g in zip(antpos.flatten(), C_phs[pi].flatten(), C_amp[pi].flatten()):
            dly[curday][pi][i] = dly[curday][pi].get(i,[]) + [tau]
            gain[curday][pi][i] = gain[curday][pi].get(i,[]) + [g]
    jday = curday
    f.close()
print dly.keys()

if opts.lst:
    if times[-1] - times[0] > 1:
        print 'Cannot plot against lst. Too many days. Reverting to days.'
        opts.lst = False
        pass
    else:
        lsts = n.array(lsts)
        if lsts[0] > lsts[-1]:
            diff = 24 - lsts[0]
            lsts = ((lsts + diff)%24) - diff

print antpos

fig = p.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
if opts.lst:
    ax1t= ax1.twiny()
    s = '-'
else:
    s = '.'

def printing_coeffs(avg,num,coeff_str=''):
    if num%8==0:
        coeff_str += '\n'
    coeff_str += '%7.3f'%avg
    return coeff_str
    
gain_str = ''
delay_str = ''
daylist_d = {}
daylist_g = {}

for day in days:
    print day
    for pi in dly[day]:
        print pi
        if not daylist_d.has_key(pi):
            daylist_d[pi] = {}
            daylist_g[pi] = {}
    #  for i in dly[pi]:
        for ik,i in enumerate(antpos.flatten()):
            if not daylist_d[pi].has_key(i): 
                daylist_d[pi][i] = []
                daylist_g[pi][i] = []
            c = colors[i/8]
            d = n.array(dly[day][pi][i])
            g = n.array(gain[day][pi][i])
            #d_avg = n.average(d)
            d_avg = n.median(d)
            g_avg = n.median(g)
            daylist_d[pi][i].append(d_avg)
            daylist_g[pi][i].append(g_avg)

            gain_str = printing_coeffs(g_avg, ik, coeff_str = gain_str)
            delay_str = printing_coeffs(d_avg, ik, coeff_str = delay_str)

    print daylist_d[pi][46]


def med_var(x):
    x = n.array(x)
    return n.mean(abs(x - n.median(x))**2)

for pi in daylist_d.keys():
    for ik,i in enumerate(antpos.flatten()):
        print med_var(daylist_d[pi][i])
print '\n\n'
for pi in daylist_d.keys():
    for ik,i in enumerate(antpos.flatten()):
        print med_var(daylist_g[pi][i])
        

print 'Variance of delays'

print 'Delay Coefficients'
#print delay_str 
print '\n \n'
print 'Gain Coefficients'
#print gain_str

