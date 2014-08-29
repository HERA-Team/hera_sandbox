#! /usr/bin/env python
"""
input a maximum number
return a set of files to be destilled.  This is the main logic behind the starting of files in the still

Selection criteria:
    files which have not made it to 6-RSYNC
    sorted by julian date

"""

from PDB import *
import optparse
import sys

o = optparse.OptionParser()
#o.add_option('-i', '--infile', dest='infile',type='string',default='',
#                help='input file')
#o.add_option('-d','--desc',dest='desc', type='string',
#                help='Description of operation')
o.add_option('--include_edges',action='store_true',
                 help="return all files necessary for stilling, not just the ones we know will make it to completion")
opts, args = o.parse_args()
cursor = pdb.db.cursor()

def check_is_edge(basefile):
    q="""select jd_hi,jd_lo from observations where basefile='{basefile}'""".format(basefile=basefile)
    cursor.execute(q)
    jds = unpack(cursor.fetchall())
    if None in jds[0]:return True
    return False


MAXNfiles = args[0]
#get the total list of files that fit selection criteria.
q="""select output,basefile from history where basefile not in (select basefile from history where operation='6-RSYNC' and exit_status=0) and
operation='1-RSYNC' and exit_status=0 order by output desc; """

cursor.execute(q)
files= unpack(cursor.fetchall())
for i,f in enumerate(files):
    #check that these aren't edge files
    basefile=f[1]
    if not opts.include_edges:
        if check_is_edge(basefile):
            #Only return files with two neighbors.
            continue
        #check also the the NEIGHBORS are not edge files
        jd_lo=pdb.get('jd_lo','observations','basefile',basefile)[0][0]
        basefile_lo = pdb.get('basefile','observations','JD',jd_lo)[0][0]
        jd_hi=pdb.get('jd_hi','observations','basefile',basefile)[0][0]
        basefile_hi = pdb.get('basefile','observations','JD',jd_hi)[0][0]
        if (check_is_edge(basefile_lo) or check_is_edge(basefile_hi)):
            continue
    if i>MAXNfiles:
        #Don't exceed the maximum number of files.
        break
    print basefile



