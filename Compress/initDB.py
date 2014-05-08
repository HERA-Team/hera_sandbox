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
pdb['files'].addcol('JD','float')
pdb['files'].addcol('basefile','string')
pdb['files'].addcol('filename','string','pk')
pdb['files'].addcol('md5','string')
pdb['files'].addcol('created_on','datetime')
pdb['files'].addcol('last_modified','datetime')
pdb['files'].addcol('host','string')

pdb.addtab('hosts')
pdb['hosts'].addcol('hostname','string','pk')
pdb['hosts'].addcol('IP','string')
pdb['hosts'].addcol('username','string')
pdb['hosts'].addcol('key_file','string')

pdb.addtab('history')
pdb['history'].addcol('input','string')
pdb['history'].addcol('output','string')
pdb['history'].addcol('operation','string')
pdb['history'].addcol('starttime','datetime')
pdb['history'].addcol('stoptime','datetime')
pdb['history'].addcol('host','string')
pdb['history'].addcol('log','string')
pdb['history'].addcol('exit_status','int')

pdb.addtab('observations')
pdb['observations'].addcol('JD','float','pk')
pdb['observations'].addcol('xx','string')
pdb['observations'].addcol('xy','string')
pdb['observations'].addcol('yx','string')
pdb['observations'].addcol('yy','string')
pdb['observations'].addcol('jd_hi','string')
pdb['observations'].addcol('jd_lo','string')
pdb['observations'].addcol('created_on','datetime')
pdb['observations'].addcol('last_modified','datetime')

pdb.addtab('orders')
pdb['orders'].addcol('filename', 'string', 'pk')
pdb['orders'].addcol('status',   'string')

if __name__ == '__main__':
    pdb.drop_tables()
    pdb.populate()
    pdb.print_schema()
