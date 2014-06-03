import unittest, random, threading, time,logging
import ddr_compress.scheduler as sch
from ddr_compress.dbi import DataBaseInterface
from ddr_compress.scheduler import FILE_PROCESSING_STAGES
import numpy as n
from sqlalchemy.orm.exc import NoResultFound
class NullAction(sch.Action):
    def _command(self): return
logger = logging.basicConfig(level=logging.INFO)
#logger = logging.getLogger().setLevel(logging.INFO)

class PopulatedDataBaseInterface(DataBaseInterface):
    def __init__(self,nobs,npols,test=True):
        DataBaseInterface.__init__(self,test=test)
        self.length = 10/60./24
        self.host = 'pot0'
        self.defaultstatus='UV_POT'
        self.Add_Fake_Observations(nobs,npols)
    def Add_Fake_Observations(self,nobs,npols):
            #form up the observation list
            obslist =[]
            jds = n.arange(0,nobs)*self.length+2456446.1234
            pols = ['xx','yy','xy','yx']
            for i,pol in enumerate(pols):
                if i>=npols:continue
                for jdi in xrange(len(jds)):
                    obslist.append({'julian_date':jds[jdi],
                                    'pol':pol,
                                    'host':self.host,
                                    'filename':'zen.{jd}.uv'.format(jd=n.round(jds[jdi],5)),
                                    'length':self.length})
                    if jdi!=0:
                        obslist[-1]['neighbor_low'] = jds[jdi-1]
                    if jdi<len(jds[:-1]):
                        obslist[-1]['neighbor_high'] = jds[jdi+1]
            obsnums = self.add_observations(obslist,status=self.defaultstatus)
    

class FakeDataBaseInterface:
    def __init__(self, nfiles=10):
        self.files = {}
        for i in xrange(nfiles):
            self.files[str(i)] = 'UV-POT'
    def get_obs_status(self, filename):
        return self.files[filename]
    def get_obs_index(self, filename): #not used
        return int(filename)
    def list_observations(self):
        files = self.files.keys()
        files.sort()
        return files
    def get_neighbors(self, filename):
        n = int(filename)
        n1,n2 = str(n-1), str(n+1)
        if not self.files.has_key(n1): n1 = None
        if not self.files.has_key(n2): n2 = None
        return (n1,n2)


#need Action's that update the db
# scheduler reads the state, decides which action to do and launches. 
#  the action sends a message to the taskserver to launch a task
# the task runs a script and then updates the db

# I need an Action that does not launch a task, it just sets the current state to Action.state
# taskserver sets the PID
        
class TestSchedulerDB(unittest.TestCase):
    def setUp(self):
        self.ntimes = 10
        self.npols = 4
        self.dbi = PopulatedDataBaseInterface(self.ntimes,self.npols,test=True)
        self.files = self.dbi.list_observations()
    def test_populated(self): #do a couple of quick checks on my db population
        obsnums =  self.dbi.list_observations()
        self.assertEqual(len(obsnums),self.ntimes*self.npols)
        self.assertEqual(len(set(obsnums)),self.ntimes*self.npols)
    def test_get_new_active_obs(self):
        s = sch.Scheduler(nstills=1, actions_per_still=1, blocksize=10)
        tic = time.time()
        s.get_new_active_obs(self.dbi)
        print "time to execute get_new_active_obs:",time.time()-tic,'s'
        self.assertEqual(len(s.active_obs),self.ntimes*self.npols)
    def test_get_action(self):
        """
        """
        obsnum = self.files[5]
        s = sch.Scheduler(nstills=1, actions_per_still=1)
        tic = time.time()
        a = s.get_action(self.dbi, obsnum, ActionClass=NullAction)
        print "time to execute get_action",time.time()-tic,'s'
        self.assertNotEqual(a, None) # everything is actionable in this test
        self.assertEqual(a.task, sch.FILE_PROCESSING_LINKS[self.dbi.defaultstatus]) # check this links to the next step        
    def test_start(self):
        self.dbi = PopulatedDataBaseInterface(3,1,test=True)
        obsnums = self.dbi.list_observations()
        class SuccessAction(sch.Action):
            def _command(me):
                me.dbi = self.dbi
                #print "Action setting {obsnum} status to {status}".format(
                #        status=me.task,obsnum=me.obs)
                me.dbi.set_obs_status(me.obs,me.task)
        def all_done():
            for obsnum in obsnums:
                try:
                    if self.dbi.get_obs_status(obsnum) != 'COMPLETE': return False
                except NoResultFound as e:
                    print f
                    print e
                    raise(NoResultFound)
        s = sch.Scheduler(nstills=1, actions_per_still=1, blocksize=10)
        t = threading.Thread(target=s.start, args=(self.dbi, SuccessAction))
        t.start()
        tstart = time.time()
        completion_time = len(FILE_PROCESSING_STAGES)*3*0.4 #0.2 s per file per step
        #print "time to completion:",completion_time,'s'
        while not all_done(): 
            if time.time() - tstart > completion_time: break
            time.sleep(1)
        s.quit()
        for obsnum in obsnums: self.assertEqual(self.dbi.get_obs_status(obsnum), 'COMPLETE')
    def test_clean_completed_actions(self):
        """
        todo
        """
        self.dbi = PopulatedDataBaseInterface(3,1,test=True)
        class SuccessAction(sch.Action):
            def _command(me):
                me.dbi = self.dbi
                me.dbi.set_obs_status(me.obs,me.task)
                print "Action has status:",me.dbi.get_obs_status(me.obs)
                return None
        s = sch.Scheduler(nstills=1, actions_per_still=1, blocksize=10)
        s.get_new_active_obs(self.dbi)
        s.update_action_queue(self.dbi, ActionClass=CompleteAction)
        a = s.pop_action_queue(0)
        s.launch_action(a)
        self.assertEqual(len(s.launched_actions[0]), 1)
        time.sleep(1)
        s.clean_completed_actions(self.dbi)
        self.assertEqual(len(s.launched_actions[0]), 0)

    def test_faulty(self):
        """
        todo
        """
        for i in xrange(1):
            dbi = FakeDataBaseInterface(10)
            class FakeAction(sch.Action):
                def __init__(self, f, task, neighbors, still):
                    sch.Action.__init__(self, f, task, neighbors, still, timeout=.01)
                def _command(self):
                    if random.random() > .5: dbi.files[self.obs] = self.task
            def all_done():
                for f in dbi.files:
                    if dbi.get_obs_status(f) != 'COMPLETE': return False
                return True
            s = sch.Scheduler(nstills=1, actions_per_still=1, blocksize=10)
            t = threading.Thread(target=s.start, args=(dbi, FakeAction))
            t.start()
            tstart = time.time()
            while not all_done() and time.time() - tstart < 20:
                #print s.launched_actions[0][0].obs, s.launched_actions[0][0].task
                #print [(a.obs, a.task) for a in s.action_queue]
                time.sleep(.1)
            s.quit()
            #for f in dbi.files:
            #    print f, dbi.files[f]
            for f in dbi.files: self.assertEqual(dbi.get_obs_status(f), 'COMPLETE')
    def test_prereqs(self):
        """
        todo
        """
        dbi = FakeDataBaseInterface(3)
        a = sch.Action('1', 'UV', ['0','2'], 0)
        self.assertTrue(a.has_prerequisites(dbi))
        for k in dbi.files: dbi.files[k] = 'CLEAN-UVC'
        a = sch.Action('1', 'ACQUIRE-NEIGHBORS', ['0','2'], 0)
        self.assertTrue(a.has_prerequisites(dbi))
        dbi.files['0'] = 'UV'
        self.assertFalse(a.has_prerequisites(dbi))
        

if __name__ == '__main__':
    unittest.main()

