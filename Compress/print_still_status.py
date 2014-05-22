#! /usr/bin/env python

"""
Prints a summary of current and past Distiller execution status.
"""

from PDB import *
import sys,optparse

o = optparse.OptionParser()
#o.add_option('--log',type=str,
#                                help='Log text. Usually saved in an env variable or catted from a file.')
opts, args = o.parse_args()
cursor = pdb.db.cursor()


print "==="*10
print "The current status of the Distiller is:"
#Summarize the hosts table
#get the pot hosts
q="""SELECT hostname FROM hosts;"""
cursor.execute(q)
hosts = unpack(cursor.fetchall())
print "   Registered Hosts"
for host in hosts:
    print host[0]
print "   "+"--"*5
q="""SELECT status FROM orders;"""
cursor.execute(q)
OrderStatus = unpack(cursor.fetchall())
print " Waiting be distilled (1-RSYNC)",len([o for o in OrderStatus if o[0]=='1-RSYNC'])
distiller_steps = map(str,range(2,6))
print " Currently being distilled",len([o for o in OrderStatus if o[0][0] in distiller_steps])
print " Waiting to be uploaded to USA (6-RSYNC)",len([o for o in OrderStatus if o[0]=='6-RSYNC'])
#get the pot hosts
q="""SELECT hostname FROM hosts where hostname like 'pot%';"""
cursor.execute(q)
pots = unpack(cursor.fetchall())
for pot in pots:
    q="""select count(*) from files where filename like '%s%%';"""%(pot[0])
    cursor.execute(q)
    potfilecount = unpack(cursor.fetchall())
    print " number of files on pot0:",potfilecount[0][0]
print "   "+"--"*5
q="select output,stoptime from history where output like 'pot0%RRE' and operation='6-RSYNC' order by stoptime limit 1;"
cursor.execute(q)
recent = unpack(cursor.fetchall())[0]
print "Last compressed file deposited onto pot:"
print recent[0],recent[1]
print "Error count in last 24 hours: ",
q="select count(*) from history where stoptime < NOW() and stoptime>DATE_SUB(NOW(),INTERVAL 1 DAY) and exit_status!=0;"
cursor.execute(q)
errcount = unpack(cursor.fetchall())[0]
print errcount[0]

"""
TODO: other useful things
average still completion time?
average transfer time?
find all int(JD)s and summarize their error count, and number of files waiting on pot0
print_still_history.py # a nice view of the history --err option only prints the erroring files and their logs


"""

