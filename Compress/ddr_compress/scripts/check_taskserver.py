#! /usr/bin/env python
import pytz
from datetime import datetime
from dateutil import parser
timeout = 15*60 #timeout in seconds
stills = ['cask0','cask1','still1','still2','still3','still4']
ucb_tz = pytz.timezone('America/Los_Angeles')
SA_tz = pytz.timezone('Africa/Johannesburg')
for host in stills:
    print host,
    logfile = '/srv/persistent/{host}/log/still/still_taskserver.log'.format(host=host)
    try:
        lines = open(logfile).readlines()
    except(IOError):
        print "[ERROR: Log not found]"
        continue
    D = None
    for line in lines:
        if line.count('alive')>0:
            D = line.split(' - ')[0].replace(',','.')
    if not D is None:
        #alivedate = datetime.strptime(D+' -0700','%Y-%m-%d %H:%M:%S,%f %z')
        alivedate = parser.parse(D+' -0700')
        if (datetime.now(SA_tz) - alivedate).total_seconds()>timeout:
            print '[FAIL]',
            print ' last log:',alivedate.astimezone(pytz.UTC).strftime('%Y-%m-%d %H:%M:%S %Z')
        else: print '[SUCCESS]'
    else:
        print "[ERROR: NO ALIVE MESSAGES]"
