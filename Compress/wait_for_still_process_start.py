#! /usr/bin/env python
"""
input: an operation (eg 2-RSYNC) and an "input" filename
Wait for an open line to exist in the history table. exit
exit status 1 for timeout
ex:
    wait_for_still_process_start.py --operation="2-RSYNC"  pot0:/data/zen.245678.4467.xx.uv
"""

from PDB import *
import optparse
import sys
from time import time,sleep
o = optparse.OptionParser()
o.add_option('--timeout',type=float,default=1.,
                        help='amount of time to wait in minutes. [default 1]')
o.add_option('--operation',type='str',
            help='the operation type we are waiting on. see http://paperwiki.berkeley.edu/doku.php/paperdb for allowable types')
opts, args = o.parse_args()
cursor = pdb.db.cursor()


starttime=time()
q="""select count(*) from history where operation='{operation}' and exit_status is NULL and output='{output}' and starttime>(NOW()-interval {timeout}
minute) and pid is NULL;""".format(
        operation=opts.operation,output=args[0],timeout=opts.timeout)
while((time()-starttime)/60.<opts.timeout):
    cursor.execute(q)
    linecount = unpack(cursor.fetchall())[0][0]
    if linecount>0:
        print "process started!"
        sys.exit()
    sleep(1)
sys.exit(1)

