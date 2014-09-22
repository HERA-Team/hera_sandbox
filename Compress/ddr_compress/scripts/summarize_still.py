#! /usr/bin/env python
"""
Prints the logs for a given obsnum or input file

"""


from ddr_compress.dbi import DataBaseInterface,gethostname,Observation,File,Log
import optparse,os,sys,re,numpy as n
import logging
from datetime import datetime,timedelta
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
    logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger('summarize_still')
dbi = DataBaseInterface()
s = dbi.Session()
print "summarizing Distiller"
OBSs = s.query(Observation)
JDs =  [OBS.julian_date for OBS in OBSs]
nights = n.sort(list(set(map(int,JDs))))
print "number of nights ingested:",len(nights)

Nobs = s.query(Observation).count()
Nprogress = s.query(Observation).filter(Observation.status!='NEW',Observation.status!='UV_POT',
        Observation.status!='NEW',Observation.status!='COMPLETE').count()
Ncomplete = s.query(Observation).filter(Observation.status=='COMPLETE').count()
print "Total observations in still:", Nobs
print "Number complete:",Ncomplete
print "Number in progress:",Nprogress
print "broken down by night [most recent activity]"
for night in nights:
    Night_complete = s.query(Observation).filter(Observation.julian_date.like(str(night)+'%'),Observation.status=='COMPLETE').count()
    Night_total = s.query(Observation).filter(Observation.julian_date.like(str(night)+'%')).count()
    OBSs = s.query(Observation).filter(Observation.julian_date.like(str(night)+'%')).all()
    obsnums = [OBS.obsnum for OBS in OBSs]
    LOG = s.query(Log).filter(Log.obsnum.in_(obsnums)).order_by(Log.timestamp.desc()).limit(1).one()
    print night,':','completeness',Night_complete,'/',Night_total,LOG.timestamp
#find all obses that have failed in the last 12 hours
FAIL_LOGs = s.query(Log).filter(Log.exit_status>0,Log.timestamp>(datetime.utcnow()-timedelta(0.5))).all()
logger.debug("found %d FAILURES"%len(FAIL_LOGs))
#break it down by stillhost
fail_obsnums = [LOG.obsnum for LOG in FAIL_LOGs]
FAIL_OBSs = s.query(Observation).filter(Observation.obsnum.in_(fail_obsnums)).all()
fail_stills = list(set([OBS.stillhost for OBS in FAIL_OBSs]))#list of stills with fails
print "fails in the last 12 hours"
for fail_still in fail_stills:
    #get failed obsnums broken down by still
    fail_count = s.query(Observation).filter(Observation.obsnum.in_(fail_obsnums),Observation.stillhost==fail_still).count()
    print fail_still,':',fail_count
sys.exit()
obsnums = [OBS.obsnum for OBS in OBSs]
still_times ={}
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
        print LOG.timestamp,
        try:print stoptime-starttime,
        except(NameError):print 'NA',
        try:
            print stoptime-computation_start
            if OBS.status == 'COMPLETE':
                try :still_times[OBS.stillhost] += [stoptime-computation_start]
                except(KeyError): still_times[OBS.stillhost] = [stoptime-computation_start]
        except(NameError):print 'NA'
print "run time summary"
print "by hosti:minutes (min,mean,max)"
for key in still_times:
    ts = [t.total_seconds() for t in still_times[key]]
    print key,':',n.round(n.min(ts)/60),n.round(n.mean(ts)/60),n.round(n.max(ts)/60)


s.close()
