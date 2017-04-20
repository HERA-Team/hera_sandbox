import os, sys, optparse, re, numpy as np
"""
Select alternating files for LST binning 
"""
o = optparse.OptionParser()
o.add_option('--odd', action='store_true', help='Toggle odd bin')
o.add_option('--verbose', action='store_true', help='Toggle verbosity.')
opts,args = o.parse_args(sys.argv[1:])

def file2jd(file):
    return float(re.findall('\D+(\d+.\d+)\D+',file)[0])
for i,filename in enumerate(args):
    if not opts.odd:
        if i%2==0:
            print filename
    else:
        if not i%2==0:
            print filename
