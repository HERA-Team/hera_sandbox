#! /usr/bin/env python
"""
Returns true if there are any processes in pdb.orders whose status isn't 1-RSYNC or 6-RSYNC (i.e. no processes are running on the still.
Returns false otherwise.

DFM
"""

from PDB import *

IdleStatus = ['1-RSYNC', '6-RSYNC']

cursor = pdb.db.cursor()
q="""SELECT status FROM orders;"""
cursor.execute(q)
OrderStatus = unpack(cursor.fetchall())

isIdle = True
for os in OrderStatus:
    # if there is something like 2-RSYNC, you should see false.
    isIdle = isIdle and (os[0] in IdleStatus)

if isIdle:
    print 'true'
else:
    print 'false'
