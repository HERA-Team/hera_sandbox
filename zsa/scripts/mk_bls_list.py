#! /usr/bin/env python
from numpy import array
import optparse,sys

o = optparse.OptionParser()
o.add_option('--fname', '-f', action='store',
            help='Name of file (with full path) to store list.')
o.add_option('--sep', '-s', action='store',
            help='Separation type given as a list. e.g. 0,1 for 0 row spacing and 1 column spacing.')
opts,args = o.parse_args(sys.argv[1:])

ant = array([[49, 41, 47, 19, 29, 28, 34, 51],
             [10,  3, 25, 48, 24, 55, 27, 57],
             [ 9, 58,  1,  4, 17, 13, 56, 59],
             [22, 61, 35, 18,  5, 32, 30, 23],
             [20, 63, 42, 37, 40, 14, 54, 50],
             [43,  2, 33,  6, 52,  7, 12, 38],
             [53, 21, 15, 16, 62, 44,  0, 26],
             [31, 45,  8, 11, 36, 60, 39, 46]])

def write_sep(name, sep=(0,1)):
    f = open(name,'w')
    st = ''
    if sep[0] < 0: start=sep[0]*-1
    else : start=0
    for i in range(ant.shape[1] - sep[1]):
        for j in range(start,ant.shape[0]-(sep[0]+start)):
            bl = str(ant[j,i]) + '_' + str(ant[j+sep[0],i+sep[1]]) + ','
            st += bl
    #leave out last comma.
    print st[:-1]
    f.write(st[:-1])

def parse_sep(s):
    '''Parses bl separation type given from options.'''
    return map(int,s.split(','))
     
write_sep(opts.fname, parse_sep(opts.sep))
