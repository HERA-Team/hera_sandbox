#!/usr/bin/env python
"""
Move all data in the specified directory into daily subdirectories. Day defined as the round JD. This is fine for SA GMT+2 and GB GMT-4
Assumes filenames of the type
<prefix>.2455900.154677.<postfix>

The files:
zen.2455899.1234.uv
zen.2455900.5745.uv

Will go to 
psa899/zen.2455899.1234.uv
psa900/zen.2455900.5745.uv

usage:
daily_move.py /path/to/data/z*uv

"""
DEBUG = False

import sys,glob as gl,os,shutil

#make relevant paths
rootdir = sys.argv[-1]
if rootdir[-1]==os.path.sep: rootdir = rootdir[:-1]
rootdir = os.path.dirname(rootdir)

prefix = 'psa'
dir_prefix = os.path.join(rootdir,prefix)


#find files
#allfiles = gl.glob(os.path.join(rootdir,))
allfiles = sys.argv[1:]
JDs = [name.split('.')[1] for name in allfiles]
tJDs = [name.split('.')[1][-3:] for name in allfiles] #tJD = truncated Julian Date
#print tJDs,[dir_prefix+tJD for tJD in tJDs]
paths = dict(zip(allfiles,[dir_prefix+tJD for tJD in tJDs]))
#print paths

unique_tJDs = set(tJDs)
mkdirs = dict(zip(unique_tJDs,[dir_prefix+unique_tJD for unique_tJD in unique_tJDs]))
print "creating directories"
#make the target directories
for tJD in unique_tJDs:
    print mkdirs[tJD],'..',
    if not os.path.exists(mkdirs[tJD]):
	os.mkdir(mkdirs[tJD])
	print "[created]"
    else: print '[exists]'
if DEBUG: "I would have moved the following %d files"%(len(allfiles))
else: print "moving %d files"%(len(allfiles))
print "eg. ",paths.keys()[0],paths[paths.keys()[0]]
for file in paths:
    print file,paths[file]
    if not DEBUG: shutil.move(file,paths[file])
