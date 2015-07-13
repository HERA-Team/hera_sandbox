#! /usr/bin/env python
"""
Subtracts a nightly average to get to zero mean visibility
uv_avg > nightly_avg.py > sub_nightly_avg.py
"""
def file2jd(file):
    return float(re.findall('\D+(\d+.\d+)\D+',file)[0])

import aipy as a, numpy as n, sys, os, optparse, pickle,re

o = optparse.OptionParser()
o.set_usage('sub_nightly_avg.py [options] *.uv')
o.add_option('--xtalk_dir',type=str,
     help='Location to look for nightly xtalk files. Default="."')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

for file in args:
    outfile = os.path.basename(file) +'x'
    print file, ' > ', outfile
    if os.path.exists(outfile):
        print "File exists. Skipping."
        continue
    jd = file2jd(os.path.basename(file))
    night = int(jd)
    
    night_avg_file_xx = opts.xtalk_dir+str(night)+'.xx.avg.pkl'
    night_avg_file_xy = opts.xtalk_dir+str(night)+'.xy.avg.pkl'
    night_avg_file_yx = opts.xtalk_dir+str(night)+'.yx.avg.pkl'
    night_avg_file_yy = opts.xtalk_dir+str(night)+'.yy.avg.pkl'
    
    NAFstr = '%s %s %s %s'%(night_avg_file_xx,night_avg_file_xy,night_avg_file_yx,night_avg_file_yy)
    print NAFstr
    
    if not os.path.exists(night_avg_file_xx) or not os.path.exists(night_avg_file_xy) or not os.path.exists(night_avg_file_yx) or not os.path.exists(night_avg_file_yy):
        print "please create nightly averages: %s"%(NAFstr)
        print "Skipping"
        continue
    uvi = a.miriad.UV(file)
    uvo = a.miriad.UV(outfile,status='new')
    uvo.init_from_uv(uvi)
    
    F_xx = open(night_avg_file_xx)
    AVG_xx = pickle.load(F_xx)
    F_xx.close()
    
    F_xy = open(night_avg_file_xy)
    AVG_xy = pickle.load(F_xy)
    F_xy.close()
    
    F_yx = open(night_avg_file_yx)
    AVG_yx = pickle.load(F_yx)
    F_yx.close()
    
    F_yy = open(night_avg_file_yy)
    AVG_yy = pickle.load(F_yy)
    F_yy.close()
    
    def mfunc(uv,preamble,data):
    	uvw,t,(i,j) = preamble
    	bl = "%d_%d"%(i,j)
    	if i!=j: 
    		if uvi['pol']==-5: data -= AVG_xx.get(bl,0)
    		elif uvi['pol']==-7: data -= AVG_xy.get(bl,0)
    		elif uvi['pol']==-8: data -= AVG_yx.get(bl,0)
    		elif uvi['pol']==-6: data -= AVG_yy.get(bl,0)
    		else:
    			print 'INVALID POL IN UVI??'
    			sys.exit(1)
    	return preamble, data
        
    uvo.pipe(uvi,mfunc=mfunc,
        append2hist='\n sub_nightly_avg.py: %s'%(NAFstr))
