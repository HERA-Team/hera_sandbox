#! /usr/bin/env python
"""
Input a list of files and insert into the db.  The files must exist and be findable on the filesystem
NB filenames must be FULL PATH. If the root is not '/' for all files it will exit

KEY NOTE: Assumes all files are contiguous.  I sort the files by jd and then match up neighboring pols as neighbors for the
   ddr algorithm

"""


from ddr_compress.dbi import DataBaseInterface,gethostname
import optparse,os,sys,re

def file2jd(zenuv):
    return re.findall(r'\d+\.\d+', zenuv)[0]
def file2pol(zenuv):
    return re.findall(r'\.(.{2})\.',zenuv)[0]
o = optparse.OptionParser()
o.set_usage('add_observations.py *.uv')
o.set_description(__doc__)
o.add_option('--length',default=10,
        help='length of the input observations in minutes [default=10]')
opts, args = o.parse_args(sys.argv[1:])

#connect to the database
dbi = DataBaseInterface()

#check that all files exist
for filename in args:
    assert(filename.startswith('/'))
    assert(os.path.exists(filename))
#now run through all the files and build the relevant information for the db
# get the pols
pols = []
for filename in args:
    pols.append(file2pol(filename))
pols = list(set(pols))#these are the pols I have to iterate over
print "found the following pols",pols
obsinfo = []
for pol in pols:
    files = [filename for filename in args if file2pol(filename)==pol]#filter off all pols but the one I'm currently working on
    files.sort()
    for i,filename in enumerate(files):
        obsinfo.append({
            'julian_date' : float(file2jd(filename)),
            'pol'     :     file2pol(filename),
            'host' :        gethostname(),
            'filename' :    filename,
            'length'  :     opts.length/60./3600 #note the db likes jd for all time units
                })
        if i!=0:
            obsinfo[-1].update({'neighbor_low':file2jd(files[i-1])})
        if i!=(len(files)-1):
            obsinfo[-1].update({'neighbor_high':file2jd(files[i+1])})
assert(len(obsinfo)==len(args))
print "adding {len} observations to the still db".format(len=len(obsinfo))
dbi.add_observations(obsinfo)
dbi.test_db()
print "done"



