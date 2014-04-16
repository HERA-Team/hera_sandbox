#! /usr/bin/env python

"""
Creates a database object for parsing and descending into the paper database.
>>> from initDB import pdb

If run as a script, this will populate and print the schema of a new database.
$ ./initDB.py

DFM
"""

import MySQLdb

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
    "float"   : "DOUBLE",
    "datetime": "TIMESTAMP"
}

class column(object):
    """
    Container for information to a single column in a database.
    """
    def __init__(self, name, dtype, key=None):
        self.name = name
        self.dtype = dtype
        self.key = key

    def parse_key(self):
        """
        For initializing pdb. Adds the proper syntax for denoting the primary key of a table.
        """
        try:
            if self.key.startswith("pk"):
                return " PRIMARY KEY,"
            else:
                return ","
        except(AttributeError):
            return ","

    def init(self):
        """
        Add an entry to the INSERT statment for the column using proper mysql syntax.
        """
        q = " %s %s"%(self.name, sqltypes[self.dtype])
        q += self.parse_key()
        return q

    def link(self,tabname):
        """
        Adds a link (foreign key) to the column.
        """
        try:
            if self.key.startswith("fk"):
                ref = self.key.split(':')[-1]
                return "ALTER TABLE %s ADD FOREIGN KEY (%s) REFERENCES %s;"%(tabname, self.name, ref)
            else:
                return None
        except(AttributeError):
            return None

class table(object):
    """
    An object for information of a table in pdb.
    """
    def __init__(self, name):
        self.name = name
        self.cols = {}

    def __getitem__(self, key):
        return self.cols[key]

    def addcol(self, name, dtype, key=None):
        """
        Adds a column to the table object (Not the db). If it's the primary key, table.pk gets populated with the column name.
        """
        self.cols[name] = column(name,dtype,key)
        if key=='pk':
            self.pk = name

    def init(self):
        """
        Write the command for creating the table in the database.
        """
        q = """CREATE TABLE IF NOT EXISTS %s ("""%self.name
        for c in self.cols:
            q += self[c].init()
        q = q[:-1] + ");"
        return q

    def link(self):
        """
        Generate mysql command for adding links among tables.
        """
        q = ""
        for c in self.cols:
            _q = self[c].link(self.name)
            if not _q is None:
                q += _q
        return q

class db(object):
    """
    An object to handle the paper database.
    """
    def __init__(self,name):
        self.name = name
        self.tabs = {}
        self.db = MySQLdb.connect(HOST, USER, PASSWD, DBNAME)
        self.db.autocommit(True)

    def __getitem__(self, key):
        return self.tabs[key]

    def __del__(self):
        self.db.close()

    def addtab(self,name):
        """
        Add a table to the database object --- this doesn't add a new table to the actual database.
        """
        self.tabs[name] = table(name)

    def print_schema(self):
        """
        Sends a human-readable schema of the database to stdout.
        """
        cout = "%s\n"%self.name
        for tname in self.tabs.keys():
            cout += " --- %s\n"%tname
            for col in self[tname].cols:
                cout += "\t --- %s (%s)"%(self[tname][col].name,self[tname][col].dtype)
                if not self[tname][col].key is None:
                    cout += " [%s]"%self[tname][col].key
                cout += "\n"
        print cout

    def initialize(self):
        """
        Populate an empty database with properly-linked tables.
        """
        #first create the tables:
        for t in self.tabs.keys():
            q = self[t].init()
            self.db.query(q)
        #next link foreign keys:
        for t in self.tabs.keys():
            q = self[t].link()
            for _q in q.split(';'):
                #a wrapper to deal with _mysql's inability to parse multiple commands on the same line.
                if not _q == "":
                    _q += ";"
                    self.db.query(_q)

    def drop_tables(self):
        """
        Deletes all tables from the database. Good for testing.
        """
        for t in self.tabs:
            q = "DROP TABLE IF EXISTS %s;"%t
            print q
            self.db.query(q)

    def has_record(self, tabname, primarykey):
        """
        Returns true if table pdb.tabname contains a row whose entry for the primary key is given by primarykey.
        """
        cursor = self.db.cursor()
        q = "SELECT EXISTS(SELECT 1 FROM %s WHERE %s='%s');"%(tabname, self[tabname].pk, primarykey)
        cursor.execute(q)
        result = cursor.fetchone()[0]
        return bool(int(result))

    def addrow(self, tabname, values):
        """
        Adds a row to pdb.tabname whose column/value pairs are given as the key/value pairs of dictionary 'values'
        """
        q = """INSERT INTO %s ("""%tabname
        q += ", ".join(values.keys())
        q += ") VALUES ("
        q += ", ".join(self.format_values(tabname, values))
        q += ");"
        print q
        self.db.query(q)

    def format_values(self, tabname, v):
        """
        Converts python strings into something that mysql understands.
        """
        vret = []
        for vi in v.keys():
            if self[tabname][vi].dtype == 'string':
                vret.append("'%s'"%v[vi])
            elif self[tabname][vi].dtype == 'float':
                vret.append(v[vi])
            elif self[tabname][vi].dtype == 'datetime':
                vret.append(v[vi])
        return vret

    def get(self, target, tab, col, val):
        """
        retrieve target column of entries in table tab, whose column is equal to val.
        """
        q = """SELECT %s FROM %s WHERE %s=%s;"""%(target, tab, col, self.format_values(tab, {col:val})[0])
        print q
        cursor = self.db.cursor()
        cursor.execute(q)
        return cursor.fetchone()

    def update(self, target_col, target_val, tab, col, val):
        """
        Change the entry of target_col to target_val in table tab for rows with col=val.
        """
        target_val = self.format_values(tab, {target_col: target_val})[0]
        val = self.format_values(tab, {col: val})[0]
        q = """UPDATE %s SET %s=%s WHERE %s=%s"""%(tab, target_col, target_val, col, val)
        self.db.query(q)

#manually enter schema

pdb = db(DBNAME)

pdb.addtab('files')
pdb['files'].addcol('JD','float')#,'fk:observations(JD)')
pdb['files'].addcol('basefile','string')#,'fk:files(basefile)')
pdb['files'].addcol('filename','string','pk')
pdb['files'].addcol('md5','string')
pdb['files'].addcol('created_on','datetime')
pdb['files'].addcol('last_modified','datetime')
pdb['files'].addcol('host','string','fk:hosts(hostname)')

pdb.addtab('hosts')
pdb['hosts'].addcol('hostname','string','pk')
pdb['hosts'].addcol('IP','string')
pdb['hosts'].addcol('username','string')
pdb['hosts'].addcol('key_file','string')

pdb.addtab('history')
pdb['history'].addcol('input','string', 'fk:files(filename)')
pdb['history'].addcol('output','string','fk:files(filename)')
pdb['history'].addcol('operation','string')
pdb['history'].addcol('timestamp','datetime')

pdb.addtab('observations')
pdb['observations'].addcol('JD','float','pk')
pdb['observations'].addcol('xx','string','fk:files(filename)')
pdb['observations'].addcol('xy','string','fk:files(filename)')
pdb['observations'].addcol('yx','string','fk:files(filename)')
pdb['observations'].addcol('yy','string','fk:files(filename)')
pdb['observations'].addcol('jd_hi','string')
pdb['observations'].addcol('jd_lo','string')
pdb['observations'].addcol('created_on','datetime')
pdb['observations'].addcol('last_modified','datetime')


if __name__ == '__main__':
    pdb.initialize()
    pdb.print_schema()
