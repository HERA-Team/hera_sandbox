#! /usr/bin/env python
"""
Creates a database object for parsing and descending into the paper database.
>>> from initDB import pdb

DFM
"""

import MySQLdb
import hashlib

#global definitions
TEST=True
HOSTIP='10.0.1.20'
USER='obs'
PASSWD='P9ls4R*@'
if TEST:
    DBNAME='test'
else:
    DBNAME='psa128'

#construct objects to lay out database schema.

def unpack(xin):
    """
    descends into nested tuples and recasts as list
    """
    if xin==None:
        return []
    else:
        xout = []
        for i in xin:
            if not type(i) == tuple:
                xout.append(i)
            else:
                xout.append(unpack(i))
        return xout

def gethostname():
    from subprocess import Popen,PIPE
    hn = Popen(['bash','-cl','hostname'], stdout=PIPE).communicate()[0].strip()
    return hn

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
        q = " %s %s"%(self.name, self.dtype)
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
    def __init__(self,name,verbose=False):
        self.name = name
        self.verbose = verbose
        self.tabs = {}
        self.db = MySQLdb.connect(HOSTIP, USER, PASSWD, DBNAME)
        self.db.autocommit(True)
        for tab in self.get_tables():
            self.addtab(tab)
            for col,dtype,key in self.get_cols(tab):
                self[tab].addcol(col,dtype,key=key)

    def __getitem__(self, key):
        return self.tabs[key]

    def __del__(self):
        self.db.close()

    def get_tables(self):
        """
        return a list of table names in pdb.
        """
        c = self.db.cursor()
        c.execute("SHOW TABLES;")
        return [t[0] for t in unpack(c.fetchall())]

    def get_cols(self, table):
        """
        return a list of column name/datatype pairs for table 'table'
        """
        c = self.db.cursor()
        c.execute("SHOW COLUMNS IN %s;"%table)
        _cols = c.fetchall()
        cols = []
        for name,dtype,null,key,default,extra in _cols:
            if key == "PRI":
                key = "pk"
            else:
                key = None
            cols.append([name,dtype,key])
        return cols

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

    def populate(self):
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
                    if self.verbose: print _q
                    self.db.query(_q)

    def drop_tables(self):
        """
        Deletes all tables from the database. Good for testing.
        """
        for t in self.tabs:
            q = "DROP TABLE IF EXISTS %s;"%t
            if self.verbose: print q
            self.db.query(q)

    def has_record(self, tabname, key, col=None):
        """
        Returns true if table pdb.tabname contains a row whose entry for the column 'col' given by key. If no column is given, the primary key is
        assumed.
        """
        cursor = self.db.cursor()
        if col is None:
            q = "SELECT EXISTS(SELECT 1 FROM %s WHERE %s='%s');"%(tabname, self[tabname].pk, key)
        else:
            q = "SELECT EXISTS(SELECT 1 FROM %s WHERE %s='%s');"%(tabname, col, key)
        if self.verbose: print q
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
        if self.verbose: print q
        self.db.query(q)

    def delrow(self, tabname, pk):
        """
        deletes a record from pdb.tabname whose primary key is pk
        """
        pk = self.format_values(tabname, {self[tabname].pk:pk})[0]
        q = """DELETE FROM %s WHERE %s=%s"""%(tabname, self[tabname].pk, pk)
        if self.verbose: print q
        self.db.query(q)

    def format_values(self, tabname, v):
        """
        Converts python strings into something that mysql understands.
          --- tabname is the table you're updating
          --- v is a dictionary with {'column name': column value} pairs.
              ---- example
              ---- >>> python_hostname = "qmaster"
              ---- >>> v = {'hostname': python_hostname}
        """
        vret = []
        for vi in v.keys():
            if self[tabname][vi].dtype in ['varchar(256)','mediumtext']:
                vret.append("'%s'"%v[vi])
            else:
                vret.append(v[vi])
        return vret

    def get(self, target, tab, col, val):
        """
        retrieve target column of entries in table tab, whose column is equal to val.
        """
        q = """SELECT %s FROM %s WHERE"""%(target,tab)
        if type(col) is list:
            constraints = [" %s=%s "%(c, self.format_values(tab, {c:v})[0]) for (c,v) in zip(col,val)]
            q += 'and'.join(constraints)[:-1]+';'
        else:
            q+=" %s=%s;"%(col, self.format_values(tab, {col:val})[0])
        if self.verbose: print q
        cursor = self.db.cursor()
        cursor.execute(q)
        return unpack(cursor.fetchall())

    def update(self, target_col, target_val, tab, col, val):
        """
        Change the entry of target_col to target_val in table tab for rows with col=val.
        """
        target_val = self.format_values(tab, {target_col: target_val})[0]
        val = self.format_values(tab, {col: val})[0]
        q = """UPDATE %s SET %s=%s WHERE %s=%s;"""%(tab, target_col, target_val, col, val)
        if self.verbose: print q
        self.db.query(q)

pdb = db(DBNAME,verbose=False)

def md5sum(fname):
    """
    calculate the md5 checksum of a file whose filename entry is fname.
    """
    fname = fname.split(':')[-1]
    BLOCKSIZE=65536
    hasher=hashlib.md5()
    try:
        afile=open(fname, 'rb')
    except(IOError):
        afile=open("%s/visdata"%fname, 'rb')
    buf=afile.read(BLOCKSIZE)
    while len(buf) >0:
        hasher.update(buf)
        buf=afile.read(BLOCKSIZE)
    return hasher.hexdigest()

if __name__ == "__main__":
    pdb.print_schema()
