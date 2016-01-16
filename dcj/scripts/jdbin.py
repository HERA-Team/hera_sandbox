import sys,optparse,re,os
import numpy as n
"""
count the number of files in a night
"""
o = optparse.OptionParser()
o.set_usage('jdbin.py [options] *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

def file2jd(file):
    return float(re.findall('\D+(\d+.\d+)\D+',file)[0])
jds = []
for filename in args:
    basename = os.path.basename(filename)
    jds.append(float(file2jd(basename)))
jds = n.array(jds)
jdfloor = n.array([n.floor(jd) for jd in jds])
days = list(set(jdfloor))
for day in days:    
    print int(day),n.sum(jdfloor==day)
