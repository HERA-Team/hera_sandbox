#! /usr/bin/env python
"""
input a host
return the number of currently executing threads on that host
"""

from PDB import *
import optparse
import sys

o = optparse.OptionParser()
#o.add_option('-i', '--infile', dest='infile',type='string',default='',
#                help='input file')
#o.add_option('-d','--desc',dest='desc', type='string',
#                help='Description of operation')
opts, args = o.parse_args()
cursor = pdb.db.cursor()
hostname=args[0]

q="""select count(*) from history where exit_status=NULL and host='{hostname}';""".format(hostname=hostname)
cursor.execute(q)
threadcount= unpack(cursor.fetchall())
print threadcount[0][0]
