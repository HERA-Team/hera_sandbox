#! /usr/bin/env python
import sys, os, math, optparse
o = optparse.OptionParser()
o.add_option('-w', '--wrap', dest='wrap', action='store_true',
    help='Instead of pulling nothing for indices off the end of the list, wrap around and repull arguments from the beginning.')
o.add_option('-n', '--num', dest='num', type='int',
    help='Number of files to group together')
opts,args = o.parse_args(sys.argv[1:])

chunk = {}
for i in xrange(len(args)/opts.num):
    chunk[i] = ''
    for j in xrange(opts.num):
        #print i,j
        chunk[i]+=args[(i*opts.num)+j]+','

for block in chunk.values():
    #print block
    os.system("combine_freqs.py -n 1024 -u {%s}" % block)
