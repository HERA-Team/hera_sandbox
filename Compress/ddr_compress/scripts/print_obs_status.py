#! /usr/bin/env python

from ddr_compress.dbi import DataBaseInterface,Observation,File
import optparse,os,sys


o = optparse.OptionParser()
o.set_usage('print_obs_status [options] zen.jd.uv')
o.set_description(__doc__)
opts, args = o.parse_args(sys.argv[1:])


dbi = DataBaseInterface()
print "thing    filename    status   stillhost"
for filename in args:
    #get the obsid of the filename
    filename = os.path.basename(filename)
    s = dbi.Session()
    OBS = s.query(Observation).filter(Observation.files.any(File.filename.like('%{filename}'.format(filename=filename)))).one()
    print 'INPUT            ',filename,OBS.status,OBS.stillhost
    print '  - High neighbor',
    try:
        high = OBS.high_neighbors[0]
        print os.path.basename(high.files[0].filename),high.status
    except(IndexError):
        print ' None'
    print '  - Low neighbor ',
    try:
        low = OBS.low_neighbors[0]
        print os.path.basename(low.files[0].filename),low.status
    except(IndexError):
        print ' None'
    s.close()
