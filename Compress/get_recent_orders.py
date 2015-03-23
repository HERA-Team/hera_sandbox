#! /usr/bin/env python
"""
Return a space-separated (i.e. BASH-readable) list of filenames who are
    1) listed in the orders table,
    2) have the most recent (integer) JD timestamp.

DFM
"""

from PDB import *
import re

cursor = pdb.db.cursor()
#first, find the most recent JD
q1="""SELECT basefile FROM orders WHERE (NOT status='6-RSYNC') and (NOT status='0-CREATED') ORDER BY basefile DESC;"""
cursor.execute(q1)
try:
    fname = unpack(cursor.fetchone())[0]
    jd_int = re.findall(r'\d+\.',fname)[0]
    #next pull all files with that JD
    #q2="SELECT %s FROM %s WHERE %s='%s' AND filename REGEXP '\%%s\%'"%('filename','orders','status',sys.argv[1],jd_int)
    files = []
    for pi in 'xy':
        for pj in 'xy':
            q2="SELECT basefile FROM orders WHERE (basefile LIKE '%{}%{}%');".format(jd_int,pi+pj)
            cursor.execute(q2)
            files +=[ c[0] for c in unpack(cursor.fetchall())]
    print ' '.join(files)
except(IndexError):
    print ""
