#! /usr/bin/env python
"""
Input a list of files and insert into the db.  The files must exist and be findable on the filesystem
NB filenames must be FULL PATH. If the root is not '/' for all files it will exit

KEY NOTE: Assumes all files are contiguous.  I sort the files by jd and then match up neighboring pols as neighbors for the
   ddr algorithm

"""


from ddr_compress.dbi import DataBaseInterface,gethostname,Observation,File
import optparse,os,sys,re,numpy as n
import logging
def file2jd(zenuv):
    return re.findall(r'\d+\.\d+', zenuv)[0]
def file2pol(zenuv):
    return re.findall(r'\.(.{2})\.',zenuv)[0]
o = optparse.OptionParser()
o.set_usage('reset_observations.py *.uv')
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
for filename in args:
    logger.debug("looking for file {filename}".format(filename=filename))
    s = dbi.Session()
    File = s.query(File).filter(File.filename==filename).one()#XXX note assumes we are not noting that this file is copied.
    obsnum = File.obsnum
    logger.debug("found obsnum {obsnum}".format(obsnum=obsnum))
    s.close()
    logger.debug("setting status to {status}".format(status=opts.status))
    dbi.set_obs_status(obsnum,opts.status)
    dbi.set_obs_pid(obsnum,None)




