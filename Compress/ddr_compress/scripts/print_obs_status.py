#! /usr/bin/env python

from ddr_compress.dbi import DataBaseInterface,Observation,File
import optparse,os,sys


o = optparse.OptionParser()
o.set_usage('print_obs_status [options] zen.jd.uv')
o.set_description(__doc__)
o.add_option('--obsnum',type=int,
        help="override the input and use this obsnum instead")
opts, args = o.parse_args(sys.argv[1:])


dbi = DataBaseInterface()
s = dbi.Session()
print "thing    filename    status   stillhost   obsnum"
if not opts.obsnum is None:
    OBS = s.query(Observation).filter(Observation.obsnum==opts.obsnum).one()
    print 'INPUT      ',os.path.basename(OBS.files[0].filename),OBS.status,OBS.stillhost,OBS.obsnum
    print '  - High neighbor',
    try:
        high = OBS.high_neighbors[0]
        print os.path.basename(high.files[0].filename),high.status,high.stillhost,high.obsnum
    except(IndexError):
        print ' None'
    print '  - Low neighbor ',
    try:
        low = OBS.low_neighbors[0]
        print os.path.basename(low.files[0].filename),low.status,low.stillhost,low.obsnum
    except(IndexError):
        print ' None'
else:
    for filename in args:
        #get the obsid of the filename
        filename = os.path.basename(filename)
        OBS = s.query(Observation).filter(Observation.files.any(File.filename.like('%{filename}'.format(filename=filename)))).one()
        print 'INPUT            ',filename,OBS.status,OBS.stillhost,OBS.obsnum
        print '  - High neighbor',
        try:
            high = OBS.high_neighbors[0]
            print os.path.basename(high.files[0].filename),high.status,high.stillhost,high.obsnum
        except(IndexError):
            print ' None'
        print '  - Low neighbor ',
        try:
            low = OBS.low_neighbors[0]
            print os.path.basename(low.files[0].filename),low.status,low.stillhost,low.obsnum
        except(IndexError):
            print ' None'
s.close()
