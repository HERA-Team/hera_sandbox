#! /usr/bin/env python
"""
Input a list of basefiles.
return true if all files are in an error or 6-RSYNC state
return false if some files are still going
"""
from PDB import *
import sys,optparse
from datetime import timedelta
o = optparse.OptionParser()
#o.add_option('--log',type=str,
#                                help='Log text. Usually saved in an env variable or catted from a file.')
opts, args = o.parse_args()
cursor = pdb.db.cursor()
basefiles = ','.join(["'{basefile}'".format(basefile=basefile) for basefile in args])
#count completeds+errors
q="""select count(*) from history where basefile in ({basefiles}) and (operation='6-RSYNC' and exit_status=0) or
exit_status!=0;""".format(basefiles=basefiles)
cursor.execute(q)
if unpack(cursor.fetchall())[0][0]==len(basefiles):print "true"
else:print "false"
