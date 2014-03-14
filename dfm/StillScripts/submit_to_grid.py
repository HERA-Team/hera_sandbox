#! /usr/bin/env python

import sys,os,optparse

o = optparse.OptionParser()
o.add_option('-t',dest='t',default=1)
o.add_option('-o',dest='o',default='/home/obs/Share/grid_output')
o.add_option('-q',dest='q',default='main')
opts,args = o.parse_args(sys.argv[1:])

#this is a hack to handle the still1 failure.
if opts.q == 'S1':
    opts.q = 'cask0'

command = 'qsub -t 1:%s -l h_vmem=8G -q %s.q -o %s -j y /home/obs/Share/redux/compress_one_file.sh'
command = command%(opts.t, opts.q, opts.o)
os.system(command)
