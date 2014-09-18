#! /usr/bin/env python
"""
Input a list of files and insert into the db.  The files must exist and be findable on the filesystem
NB filenames must be FULL PATH. If the root is not '/' for all files it will exit

KEY NOTE: Assumes all files are contiguous.  I sort the files by jd and then match up neighboring pols as neighbors for the
   ddr algorithm

"""


from ddr_compress.dbi import DataBaseInterface,gethostname,jdpol2obsnum
import optparse,os,sys,re,numpy as n

def file2jd(zenuv):
    return re.findall(r'\d+\.\d+', zenuv)[0]
def file2pol(zenuv):
    return re.findall(r'\.(.{2})\.',zenuv)[0]
o = optparse.OptionParser()
o.set_usage('add_observations.py *.uv')
o.set_description(__doc__)
o.add_option('--length',type=float,
        help='length of the input observations in minutes [default=average difference between filenames]')
o.add_option('-t',action='store_true',
       help='Test. Only print, do not touch db')
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
jds = []
for filename in args:
    pols.append(file2pol(filename))
    jds.append(float(file2jd(filename)))
jds = n.array(jds)
if not opts.length is None:
    djd =  opts.length/60./24
else:
    jds_onepol = n.sort([jd for i,jd in enumerate(jds) if pols[i]==pols[0]])
    djd = n.mean(n.diff(jds_onepol))
    print "setting length to ",djd,' days'
pols = list(set(pols))#these are the pols I have to iterate over
print "found the following pols",pols
obsinfo = []
for pol in pols:
    files = [filename for filename in args if file2pol(filename)==pol]#filter off all pols but the one I'm currently working on
    files.sort()
    for i,filename in enumerate(files):
        try:
            dbi.get_obs(jdpol2obsnum(float(file2jd(filename)),file2pol(filename),djd))
            print filename, "found in db, skipping"
        except:
            obsinfo.append({
                'julian_date' : float(file2jd(filename)),
                'pol'     :     file2pol(filename),
                'host' :        gethostname(),
                'filename' :    filename,
                'length'  :     djd #note the db likes jd for all time units
                    })
for i,obs in enumerate(obsinfo):
    filename = obs['filename']
    if i!=0:
        if n.abs(obsinfo[i-1]['julian_date']-obs['julian_date'])/djd<0.01:
            obsinfo[i].update({'neighbor_low':file2jd(files[i-1])})
    if i!=(len(obsinfo)-1):
        if n.abs(obsinfo[i+1]['julian_date']-obs['julian_date'])/djd<0.01:
            obsinfo[i].update({'neighbor_high':file2jd(files[i+1])})
#assert(len(obsinfo)==len(args))
if opts.t:
    print "NOT ADDING OBSERVATIONS TO DB"
    print "HERE is what would have been added"
    for obs in obsinfo:
        print obs['filename'],jdpol2obsnum(obs['julian_date'],obs['pol'],obs['length'])
elif len(obsinfo)<0:
    print "adding {len} observations to the still db".format(len=len(obsinfo))
    dbi.add_observations(obsinfo)
    dbi.test_db()
print "done"



