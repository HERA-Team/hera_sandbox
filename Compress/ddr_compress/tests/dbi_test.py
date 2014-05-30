import unittest, random, threading, time
import ddr_compress.scheduler as sch
from ddr_compress.dbi import Base,File,Observation
from ddr_compress.dbi import databaseinterface,jdpol2obsnum,obsnum2jdpol
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
import numpy as n

class TestDBI(unittest.TestCase):
    def setUp(self):
        """
        create an in memory DB and open a connection to it
        """
        self.dbi = databaseinterface(test=True)
        self.session = self.dbi.Session()
        self.jd = 2456785.123456
        self.pol = 'xx'
        self.filename='/data0/zen.2456785.123456.uv'
        self.host = 'pot0'
    def test_obsnum(self):
        obsnum = jdpol2obsnum(self.jd,self.pol)
        outjd,outpol = obsnum2jdpol(obsnum)
        self.assertEqual(self.pol,outpol)
        self.assertEqual(self.jd,outjd)

    def test_Observation_and_File(self):
        """
        Create an observation record
        """
        obsnum = jdpol2obsnum(self.jd,self.pol)
        self.observation1 = Observation(
                    julian_date=self.jd,
                    pol=self.pol,
                    status='UV-POT')
        self.file1 = File(
                    filename=self.filename,
                    host=self.host)
        self.file1.observation=self.observation1
        self.session.add(self.file1)
        self.session.add(self.observation1)
        self.session.commit()
        
        OBS = self.session.query(Observation).one()
        #check that we got our observation back
        self.assertEqual(OBS.obsnum,self.observation1.obsnum)
        #check that the obs is linked to the file
        self.assertEqual(OBS.files[0].host,self.file1.host)
        self.assertEqual(OBS.obsnum,obsnum)
    def test_add_observation(self):
        """
        use the dbi to create a record.
        basically tests the same as test_Observation_and_file
        but with the dbi wrapper
        """
        obsnum = self.dbi.add_observation(self.jd,self.pol)
        OBS = self.session.query(Observation).filter(Observation.obsnum==obsnum).one()
        self.assertEqual(OBS.julian_date,self.jd)
    def test_add_file(self):
        """
        """
        #first add the observation
        obsnum = self.dbi.add_observation(self.jd,self.pol)
        #then add a file to it
        filenum = self.dbi.add_file(obsnum,self.host,self.filename)
        #then grab the file record back
        FILE = self.session.query(File).filter(File.filenum==filenum).one()
        #and check that its right
        self.assertEqual(FILE.filename,self.filename)
    def test_add_files(self):
        #first add the observation
        obsnum = self.dbi.add_observation(self.jd,self.pol)

        files = ['/data0/zen.2456785.123456.uv','/data0/zen.2456785.123456.uvc','/data0/zen.2456785.323456.uvcR']
        for filename in files:
            filenum = self.dbi.add_file(obsnum,self.host,filename)
        #how I get all the files for a given obs
        OBS = self.session.query(Observation).filter(Observation.obsnum==obsnum).one()
        self.assertEqual(len(OBS.files),len(files))#check that we got three files
    def test_set_obs_status(self):
        """
        set the status with the dbi function then check it with 
        under the hood stuff
        """
        #first create an observation in the first place
        obsnum = self.dbi.add_observation(self.jd,self.pol)
        # then set the status to something else
        self.dbi.set_obs_status(obsnum,'UV')
        # get the status back out 
        OBS = self.session.query(Observation).filter(Observation.obsnum==obsnum).one()
        self.assertEqual(OBS.status,'UV')
    def test_get_obs_status(self):
        """
        set the status with the dbi function (cause I tested it with more basic tests already)
        """
        #first create an observation in the first place
        obsnum = self.dbi.add_observation(self.jd,self.pol)
        # then set the status to something else
        self.dbi.set_obs_status(obsnum,'UV')
        #then get the status back
        status = self.dbi.get_obs_status(obsnum)
        self.assertEqual(status,'UV')

    def test_get_neighbors(self):
        """
        First set up a likely triplet of observations
        """
        dt = 10/60./24  
        jds = n.arange(0,3)*dt+2456446.1234
        obsnums=[]
        for jd in jds:
            obsnum = self.dbi.add_observation(jd,self.pol,length=dt)
            obsnums.append(obsnum)
        neighbors = self.dbi.get_neighbors(obsnums[1])
        self.assertEqual(len(neighbors),2)
        self.assertEqual(neighbors[0],obsnums[0])
        self.assertEqual(neighbors[1],obsnums[2])
        
        

if __name__=='__main__':
    unittest.main()

