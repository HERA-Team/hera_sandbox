#
#  CAPOdB.py
#  
#
#  Created by Danny Jacobs on 3/10/09.
# 
#
"""
Module enabling databasing of calibration values and sky model.
DCJ 10 March 2009
"""
import aipy as a, numpy as n, sqlite3 as sql


class CAPO_dB():
    def __init__(self,filename='paper_test1.db'):
        self.conn = sql.connect(filename)
        self.conn.row_factory = sql.Row
        self.c = self.conn.cursor()
        try: self.c.execute('''create table aa (ant int, bp_r blob, amp  real, 
                             phsoff text, x real, y real, z real, bp_i text,
                             dly real,off real)''')
        except: pass #the db already exists or something I don't know what is wrong.
        self.c.close()
        self.aa_db_params = ['ant','bp_r','bp_i','amp','phsoff','x','y','z']
class AntennaArray(a.fit.AntennaArray,CAPO_dB):
    def __init__(self,location,ants,**kwargs):
        if kwargs.has_key('filename'): CAPO_dB.__init__(self,filename=kwargs['filename'])
        else: CAPO_dB.__init__(self)
        a.fit.AntennaArray.__init__(self, location,ants,**kwargs)
        #initialize antennae using existing table entries
        if kwargs.has_key('beam'): beam = kwargs['beam']
        if ants==[]:
           self.c.execute("""
           SELECT ant,x,y,z,dly,off,amp FROM aa 
           where rowid=(select max(rowid) from aa as a where a.ant=aa.ant)
           """)
           self.ants = []
           for db_row in self.c:
              self.ants.append(Antenna(0,0,0,beam,db_row=db_row,
              db_desc=self.c.description))
    def clear(self):
        """
        Delete _all_ entries in a table.  Use with trepidation.
        """
        c = self.conn.cursor()
        c.execute("DELETE FROM aa")
        c.close()
                        
    def save(self):
        """
        Saves parameters returned by get_params()
        """
        self.c = self.conn.cursor()
        for a in self: a.save(self.c,self.aa_db_params)
        self.conn.commit()
        self.c.close()
            
    def prm_update(self):
        """
        Gets latest antennae from db. Completely wipes and repopulates
        parameters.
        """
        self.c = self.conn.cursor()
        for a in self: a.prm_update(self.c,self.aa_db_params)
        self.c.close()
        
class Antenna(a.fit.Antenna):
    """
    Adds saving of antenna records to dB to fit.Antenna class. 
    """
    def __init__(self, x, y, z, beam, phsoff=[0.,0.], bp_r=n.array([1]),
            bp_i=n.array([0]), amp=1, pointing=(0.,n.pi/2,0), **kwargs):
        """
        Creates an antenna given a row from a db (ie db_row=cursor.fetchone()) 
        and the description field (ie db_desc=conn.description)
        """
        if kwargs.has_key('db_row'):
           db_row=kwargs['db_row']
           db_desc=kwargs['db_desc']
           #parse out the dumb description
           db_desc = [L[0] for L in db_desc]
           for pname in db_desc:
              if type(db_row[pname])==str: exec("%s = eval(%s)"%(pname,db_row[pname]))
              else: exec("%s = db_row[pname]"%(pname))
        a.fit.Antenna.__init__(self, x,y,z, beam, phsoff=phsoff,bp_r=bp_r,
            bp_i=bp_i,amp=amp)
        if kwargs.has_key('id'):self.id = kwargs['id']
        elif 'ant' in db_desc: self.id = ant
    def save(self,c,plist):
        """
        Enter the antenna's parameters as a new row in the db specified by the parent
        AntennaArray.
        save(cursor)#where cursor is a sqlite cursor.
        plist = list of parameters to save, all must have columns in db
        """
        prms = self.get_params(['*'])
        prms.update({'ant':self.id})
        names = plist
        vals = [prms[k] for k in names]
        vals = map(str,vals)
        record = "insert into aa %s values %s"%(str(tuple(names)),
                 str(tuple(vals)))
        c.execute(record)
    def prm_update(self,c,plist):
        """
        Load the antenna's most recent parameters from the db
        """
        names = str(tuple(plist))
        names = names.replace("'","")[1:-1]            
        record = """select %s from aa where ant = ? and
        rowid=(select max(rowid) from aa as a where a.ant=aa.ant)"""%(names)
        #print record
        c.execute(record,(self.id,))
        A = c.fetchone()
        dnames = [k[0] for k in c.description]
        prms = {}
        for pname in plist:
            if type(A[pname])==unicode: val = eval(A[pname].replace("u",""))
            else: val = A[pname]
           # print A[pname],type(A[pname])
            prms.update({pname:val})
        self.set_params(prms)
