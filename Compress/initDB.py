#! /usr/bin/env python

"""
If run as a script, this will populate and print the schema of a new database
$ ./initDB.py

DFM
"""

from PDB import *
#manually enter schema

pdb = db(DBNAME)

pdb.addtab('files')
pdb['files'].addcol('JD','double')
pdb['files'].addcol('basefile','varchar(256)')
pdb['files'].addcol('filename','varchar(256)','pk')
pdb['files'].addcol('md5','varchar(256)')
pdb['files'].addcol('created_on','timestamp')
pdb['files'].addcol('last_modified','timestamp')
pdb['files'].addcol('host','varchar(256)')

pdb.addtab('hosts')
pdb['hosts'].addcol('hostname','varchar(256)','pk')
pdb['hosts'].addcol('IP','varchar(256)')
pdb['hosts'].addcol('username','varchar(256)')
pdb['hosts'].addcol('key_file','varchar(256)')

pdb.addtab('history')
pdb['history'].addcol('input','varchar(256)')
pdb['history'].addcol('output','varchar(256)')
pdb['history'].addcol('operation','varchar(256)')
pdb['history'].addcol('starttime','timestamp')
pdb['history'].addcol('stoptime','timestamp')
pdb['history'].addcol('host','varchar(256)')
pdb['history'].addcol('log', 'mediumtext')
pdb['history'].addcol('exit_status','tinyint')
pdb['history'].addcol('basefile','varchar(256)')
pdb['history'].addcol('pid','int')

pdb.addtab('observations')
pdb['observations'].addcol('JD','double')
pdb['observations'].addcol('pol','varchar(256)')
pdb['observations'].addcol('basefile','varchar(256)','pk')
pdb['observations'].addcol('jd_hi','varchar(256)')
pdb['observations'].addcol('jd_lo','varchar(256)')
pdb['observations'].addcol('created_on','timestamp')
pdb['observations'].addcol('last_modified','timestamp')

pdb.addtab('orders')
pdb['orders'].addcol('basefile', 'varchar(256)', 'pk')
pdb['orders'].addcol('status',   'varchar(256)')

if __name__ == '__main__':
    pdb.drop_tables()
    pdb.populate()
    pdb.print_schema()
