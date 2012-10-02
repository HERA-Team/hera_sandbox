#! /usr/bin/env python
import sys, os, math, optparse
o = optparse.OptionParser()
o.add_option('-w', '--wrap', dest='wrap', action='store_true',
    help='Instead of pulling nothing for indices off the end of the list, wrap around and repull arguments from the beginning.')
o.add_option('-t', type='str',
    help='Same as -t qsub. ex -t 1:10 Use if you want to multitask without grid engine. Must specifiy task_id as well.')
o.add_option('--taskid',type=int,default=1,
    help='current taskid. Use in conjunction with -t for non-gridengine list splitting mode.')
o.add_option('-f',action='store_true',
    help='look in the input files for args instead [default = false]')
opts,args = o.parse_args(sys.argv[1:])

if opts.f:
    fileargs = []
    for file in args:
        lines = open(file).readlines()
        for line in lines:
            fileargs.append(line)
    args  = fileargs
if not opts.t is None:
    n = int(opts.t.split(':')[0])
    m = int(opts.t.split(':')[1])
    i = opts.taskid-1
    if (m-n) <= len(args) or not opts.wrap:
        num = int(math.ceil(float(len(args)) / (m - n + 1)))
        print ' '.join(args[num*i:num*(i+1)])
    else:
        print args[i % len(args)]
    sys.exit(0)    
try:
    n = int(os.environ['SGE_TASK_FIRST'])
    m = int(os.environ['SGE_TASK_LAST'])
    i = int(os.environ['SGE_TASK_ID']) - 1
    if (m-n) <= len(args) or not opts.wrap:
        num = int(math.ceil(float(len(args)) / (m - n + 1)))
        print ' '.join(args[num*i:num*(i+1)])
    else:
        print args[i % len(args)]
except(KeyError,ValueError): print ' '.join(args)
