#! /usr/bin/env python
"""
Prints the logs for a given obsnum or input file

"""


from ddr_compress.dbi import DataBaseInterface,gethostname,Observation,File,Log
import optparse,os,sys,re,numpy as n
import logging
def file2jd(zenuv):
    return re.findall(r'\d+\.\d+', zenuv)[0]
def file2pol(zenuv):
    return re.findall(r'\.(.{2})\.',zenuv)[0]
o = optparse.OptionParser()
o.set_usage('summarize_night.py obsnum')
o.set_description(__doc__)
#o.add_option('--length',type=float,
#        help='length of the input observations in minutes [default=average difference between filenames]')
o.add_option('-v',action='store_true',
        help='set log level to debug')
o.add_option('--status',default='UV_POT',
        help='set the observation to this status [default=UV_POT]')
opts, args = o.parse_args(sys.argv[1:])
#connect to the database
if opts.v:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('reset_observations')
dbi = DataBaseInterface()
s = dbi.Session()

OBSs = s.query(Observation).filter(Observation.julian_date.between(float(args[0]),float(args[0])+0.9999)).all()
obsnums = [OBS.obsnum for OBS in OBSs]
for i,OBS in enumerate(OBSs):
    obsnum = OBS.obsnum
    if OBS.stillhost is None:continue
    print '\t'.join(map(str,[obsnum,OBS.stillhost,dbi.get_input_file(obsnum)[2],OBS.status])),
    LOGs = s.query(Log).filter(Log.obsnum==obsnum).order_by(Log.timestamp.desc()).all()
    if len(LOGs)==0:print 'NO LOGS';continue
    stoptime=LOGs[0].timestamp
    computation_start=0
    if opts.v:
        for LOG in LOGs[1:]:
            if LOG.stage=='UV_POT' or LOG.stage=='NEW':break
            print LOG.stage+"="+str(n.round((stoptime - LOG.timestamp).total_seconds()/3600,1)),
            if LOG.exit_status != 0: print '!',
            print '({stat},{pid})'.format(stat=LOG.exit_status,pid=OBS.currentpid),
            stoptime = LOG.timestamp
        print

    else:
        for LOG in LOGs:
            logger.debug(LOG.stage+':'+str(LOG.timestamp))
            if LOG.stage=='NEW':starttime=LOG.timestamp;break
            if LOG.stage=='UV_POT': starttime=LOG.timestamp;break
            if LOG.stage=='UV': computation_start = LOG.timestamp
        try:print stoptime-starttime,
        except(NameError):print 'NA',
        try: print stoptime-computation_start
        except(NameError):print 'NA'


s.close()
