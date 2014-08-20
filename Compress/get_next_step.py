#! /usr/bin/env python
"""
This is the MAIN LOGIC Center of the distiller.
qdaemon runs this script to figure out what to do next for any file.

input: basefile.
output: taskname.  tasks it might give on output
NULL     Don't do anything
2-RSYNC  Transfer from pot? to still?
3-XRFI   correct and first-round xrfi (uvcR)
4-XRFI   invert compression, generate rfi npz file, apply to data (uvcRR, uvcRE.npz)
5-DDR    delay-delay-rate filter (uvcRRE)
6-RSYNC  Transfer from still? back to pot?
7-RSYNC  Transfer to USA
TIMEOUT  kill this current task
FINISH   This file made it all the way through and is done forever. Consider skipping.
ERROR    Something has gone horribly wrong but no processes need to be killed. Consider going back a step or possibly restarting from the beginning.

#TODO
more verbose printouts for the ERROR conditions with a -v switch.
"""

from PDB import *
import sys,optparse
from datetime import timedelta
o = optparse.OptionParser()
#o.add_option('--log',type=str,
#                                help='Log text. Usually saved in an env variable or catted from a file.')
o.add_option('--timeout',type=float,default=100,
                help='timeout in hours [default=100hrs]')
opts, args = o.parse_args()
cursor = pdb.db.cursor()

task_assignment = {
"1-RSYNC":"2-RSYNC",
"2-RSYNC":"3-XRFI",
"3-XRFI": "4-XRFI",
"4-XRFI": "5-DDR",
"5-DDR":  "6-RSYNC",
"6-RSYNC":"7-RSYNC",
"7-RSYNC":"NULL",
"RESET":"2-RSYNC"}
MAXTHREADS=8
basefile=sys.argv[1]
#first question.  Is the file currently doing anything?
q="""select operation,exit_status,starttime,stoptime from history where basefile='{basefile}' order by operation desc limit 1;""".format(basefile=basefile)
cursor.execute(q)
most_recent_task = unpack(cursor.fetchall())
if most_recent_task[0][1] is None: #if the most recent task is under way
    #impliment timeout here. not sure of units of stoptime when read out this way
    if (most_recent_task[0][3]-most_recent_task[0][2])>timedelta(hours=opts.timeout):
        print "TIMEOUT"
    else:
        print "NULL"
    sys.exit(0)

#second question.  what is the most recently completed task?
q="""select operation,host from history where stoptime is not NULL and
basefile='{basefile}' and exit_status=0 order by operation desc limit 1;""".format(basefile=basefile)
cursor.execute(q)
most_recently_completed_task,most_recent_host = unpack(cursor.fetchall())[0]
#for now just blindly give the next possible task.
#TODO: add necessary logic for each task
#for tasks on the still nodes, check the number of available threads
#for tasks needing neighbors, check that the neighbors are available
#for tasks needing an scp thread. make sure that there less than Nscps running
try:
    TASK = task_assignment[most_recently_completed_task]
except(KeyError):
    print "ERROR"
    sys.exit(0)
#BLOCKING.  Check that the system and any neighboring data sets are ready
#currently blocking scp and still executions in the same way.
#ie Available SCP threads = Available host cores
if int(TASK[0])>0:#block on the number of threads available on the host where the data are
    q="""select count(*) from history where exit_status=NULL and host='{hostname}';""".format(hostname=most_recent_host)
    cursor.execute(q)#get the number of threads on the host where this operation is supposed to happen
    threadcount= unpack(cursor.fetchall())[0][0]
    if threadcount>=MAXTHREADS:
        print "NULL"
        sys.exit(0)
if TASK in ['4-XRFI','5-DDR']:#these three tasks require neighbors. if no neighbors. wait
    #todo. Check to see if the neighbors are ordered.
    #      for now assumes they've been properly submitted by qdaemon.
    #first get the neighbor basefiles.
    jd_lo=pdb.get('jd_lo','observations','basefile',basefile)[0][0]
    basefile_lo = pdb.get('basefile','observations','JD',jd_lo)[0][0]
    jd_hi=pdb.get('jd_hi','observations','basefile',basefile)[0][0]
    basefile_hi = pdb.get('basefile','observations','JD',jd_hi)[0][0]

    #get neighboring files on the current execution host
    q="""select outfile from history where exit_status=0 and (basefile='{basefile_hi}' or basefile='{basefile_lo}') and host='{most_recent_host}' and operation='{most_recent_task}';""".format(
        basefile_hi=basefile_hi,
        basefile_lo=basefile_lo,
        most_recent_host=most_recent_host,
        most_recent_task=most_recent_task)
    cursor.execute(q)
    neighbors = unpack(cursor.fetchall())
    if len(neighbors)!=2:
        print "NULL"
        sys.exit(0)
#add any more blocks or actions HERE!
print TASK

