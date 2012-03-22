#! /usr/bin/env python
import sys, optparse
o = optparse.OptionParser()
o.add_option('--baseport', dest='baseport', type='int', default=53000,
    help='Number to add to daemon id to get port number')
o.add_option('--divstr', dest='divstr', default=',',
    help='Divider string between host:port entries.  Default ","')
opts,args = o.parse_args(sys.argv[1:])

def parseline(L):
    L = L.split()
    node = L[-3].split('@')[-1]
    port = int(L[-1])
    return '%s:%d' % (node, port+opts.baseport)
print opts.divstr.join([parseline(L) for L in sys.stdin.readlines() 
    if L.split()[0] == args[0]])
