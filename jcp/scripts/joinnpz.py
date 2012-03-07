#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n, sys, optparse

o = optparse.OptionParser()
o.set_usage('joinnpz.py file1 file2')
o.add_option('-p','--plot',dest='plot',action='store_true',
    help='Plot the files you joined. Default is False.')
             
opts,args = o.parse_args(sys.argv[1:])

file1 = n.load(args[0])
file2 = n.load(args[1])

times = n.append(file1['times'],file2['times'])
spec = n.append(file1['spec'],file2['spec'],axis=0)
afreqs = file1['afreqs']

if opts.plot:
    #p.imshow(spec.real)
    p.plot(times)
    p.show()

src,ftimes = args[0].split('__')
time_i = ftimes.split('_')[0]
time_f = args[1].split('__')[1].split('_')[1]
outname = src+'__'+time_i+'_'+time_f
print "Creating "+outname
n.savez(outname,times=times,afreqs=afreqs,spec=spec)
