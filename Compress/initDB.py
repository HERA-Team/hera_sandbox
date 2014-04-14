#! /usr/bin/env python

import _mysql

#global definitions
TEST=True
HOST='localhost'
USER='obs'
PASSWD='P9ls4R*@'
if TEST:
    DBNAME='test'
else:
    DBNAME='psa128'

#construct objects to lay out database schema.

sqltypes = {
    "string"  : "VARCHAR(256)",
    "float"   : "FLOAT",
    "datetime": "TIMESTAMP"
}

class column(object):

    def __init__(self, name, dtype, key=None):
        self.name = name
        self.dtype = dtype
        self.key = key

    def parse_key(self):
        try:
            if self.key.startswith("pk"):
                return " PRIMARY KEY,"
            else:
                return ","
        except(AttributeError):
            return ","

    def init(self):
        q = " %s %s"%(self.name, sqltypes[self.dtype])
        q += self.parse_key()
        return q

    def link(self,tabname):
        try:
            if self.key.startswith("fk"):
                ref = self.key.split(':')[-1]
                return "ALTER TABLE %s ADD FOREIGN KEY (%s) REFERENCES %s;"%(tabname, self.name, ref)
            else:
                return None
        except(AttributeError):
            return None

class table(object):

    def __init__(self, name):
        self.name = name
        self.cols = []

    def add(self, name, dtype, key=None):
        self.cols.append(column(name,dtype,key))

    def init(self):
        q = """CREATE TABLE IF NOT EXISTS %s ("""%self.name
        for c in self.cols:
            q += c.init()
        q = q[:-1] + ");"
        return q

    def link(self):
        q = ""
        for c in self.cols:
            _q = c.link(self.name)
            if not _q is None:
                q += _q
        return q

class db(object):

    def __init__(self,name):
        self.name = name
        self.tabs = {}

    def __getitem__(self, key):
        return self.tabs[key]

    def add(self,name):
        self.tabs[name] = table(name)

    def print_schema(self):
        cout = "%s\n"%self.name
        for tname in self.tabs.keys():
            cout += " --- %s\n"%tname
            for col in self[tname].cols:
                cout += "\t --- %s (%s)"%(col.name,col.dtype)
                if not col.key is None:
                    cout += " [%s]"%col.key
                cout += "\n"
        print cout

    def initialize(self, test=False):
        _db = _mysql.connect(HOST, USER, PASSWD, DBNAME)
        _db.autocommit(True)
        #first create the tables:
        for t in self.tabs.keys():
            q = self[t].init()
            print q
            _db.query(q)
        #next link foreign keys:
        for t in self.tabs.keys():
            q = self[t].link()
            for _q in q.split(';'):
                #a wrapper to deal with _mysql's inability to parse multiple commands on the same line.
                if not _q == "":
                    _q += ";"
                    print _q
                    _db.query(_q)

    def drop_tables(self):
        _db = _mysql.connect(HOST, USER, PASSWD, DBNAME)
        _db.autocommit(True)
        for t in self.tabs:
            q = "DROP TABLE IF EXISTS %s;"%t
            print q
            _db.query(q)

#manually enter schema

pdb = db(DBNAME)

pdb.add('files')
pdb['files'].add('JD','float','fk:observations(JD)')
pdb['files'].add('basefile','string')#,'fk:files(basefile)')
pdb['files'].add('filename','string','pk')
pdb['files'].add('md5','string')
pdb['files'].add('created_on','datetime')
pdb['files'].add('last_modified','datetime')
pdb['files'].add('host','string','fk:hosts(hostname)')

pdb.add('hosts')
pdb['hosts'].add('hostname','string','pk')
pdb['hosts'].add('IP','string')
pdb['hosts'].add('username','string')
pdb['hosts'].add('key_file','string')

pdb.add('history')
pdb['history'].add('input','string', 'fk:files(filename)')
pdb['history'].add('output','string','fk:files(filename)')
pdb['history'].add('operation','string')
pdb['history'].add('timestamp','datetime')

pdb.add('observations')
pdb['observations'].add('JD','float','pk')
pdb['observations'].add('xx','string','fk:files(filename)')
pdb['observations'].add('xy','string','fk:files(filename)')
pdb['observations'].add('yx','string','fk:files(filename)')
pdb['observations'].add('yy','string','fk:files(filename)')
pdb['observations'].add('jd_hi','string')
pdb['observations'].add('jd_lo','string')
pdb['observations'].add('created_on','datetime')
pdb['observations'].add('last_modified','datetime')


if __name__ == '__main__':
    pdb.initialize(test=TEST)
    pdb.print_schema()
