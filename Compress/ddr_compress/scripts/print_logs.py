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
o.set_usage('print_logs.py obsnum')
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
    logging.basicConfig(lebel=logging.INFO)
logger = logging.getLogger('reset_observations')
dbi = DataBaseInterface()

# for each file get the obsnum, then reset the status to UV_POT
obsnums = []
for arg in args:
    if arg[0]=='/': #we have a file!
        logger.debug("looking for file {filename}".format(filename=arg))
        s = dbi.Session()
        File = s.query(File).filter(File.filename==arg).one()#XXX note assumes we are not noting that this file is copied.
        obsnum = File.obsnum
        logger.debug("found obsnum {obsnum}".format(obsnum=obsnum))
        s.close()
        obsnums.append(obsnum)
    else:
        obsnums.append(arg)#it must be an obsnum
logger.debug("found %d obsnums"%len(obsnums))

#get the logs
s = dbi.Session()
#things I can print. 
#  logtext
#  exit_status
#  timestamp
#  stage
#  obsnum
for obsnum in obsnums:
    LOGs = s.query(Log).filter(Log.obsnum==obsnum).order_by(Log.timestamp.asc()).all()
    print "="*80
    filename = dbi.get_input_file(obsnum)
    filename = filename[0]+':'+filename[1]+'/'+filename[2]
    print "==    History for obsnum={obsnum} file={filename}     ==".format(
                                    obsnum = obsnum,
                                    filename=filename)
    print "=="
    print "=="
    print "=="
    for LOG in LOGs:
        print "-"*80
        print "  obsnum:{obsnum} ".format(obsnum=obsnum)
        print "  processing stage:{stage}".format(stage=LOG.stage)
        print "  timestamp: {datetime}".format(datetime=LOG.timestamp)
        if LOG.exit_status==0:
            print "  exit_status: SUCCESS"
        else:
            print "  exit_status:  FAIL [{exit_status}]".format(exit_status=LOG.exit_status)
        print 
        print 
        print LOG.logtext

s.close()

 




