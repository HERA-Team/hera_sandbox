import unittest, random, threading, time
import ddr_compress.scheduler as sch
from ddr_compress.dbi import DataBaseInterface
import numpy as n
class NullAction(sch.Action):
    def _command(self): return


class PopulatedDataBaseInterface(DataBaseInterface):
    def __init__(self,nobs,npols,test=True):
        DataBaseInterface.__init__(self,test=test)
        self.length = 10/60./24
        self.host = 'pot0'
        print "making fake db" 
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
            obsnums = self.add_observations(obslist)
    

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

class TestAction(unittest.TestCase):
    def setUp(self):
        self.files = ['1','2','3']
        self.still = 0
        self.task = 'UVC'
    def test_attributes(self):
        a = sch.Action(self.files[1], self.task, [self.files[0],self.files[2]], self.still)
        self.assertEqual(a.task, self.task)
        # XXX could do more here
    def test_priority(self):
        a = sch.Action(self.files[1], self.task, [self.files[0],self.files[2]], self.still)
        self.assertEqual(a.priority, 0)
        a.set_priority(5)
        self.assertEqual(a.priority, 5)
    def test_prereqs(self):
        a = sch.Action(self.files[1], self.task, [self.files[0],self.files[2]], self.still)
        self.assertTrue(a.has_prerequisites(None))
        # XXX more here
    def test_timeout(self):
        a = NullAction(self.files[1], self.task, [self.files[0],self.files[2]], self.still, timeout=100)
        self.assertRaises(AssertionError, a.timed_out)
        t0 = 1000
        a.launch(launch_time=t0)
        self.assertFalse(a.timed_out(curtime=t0))
        self.assertTrue(a.timed_out(curtime=t0+110))
    def test_action_cmp(self):
        priorities = range(10)
        actions = [sch.Action(self.files[1], self.task, [self.files[0],self.files[2]], self.still) for p in priorities]
        random.shuffle(priorities)
        for a,p in zip(actions,priorities): a.set_priority(p)
        actions.sort(cmp=sch.action_cmp)
        for cnt,a in enumerate(actions):
            self.assertEqual(a.priority, cnt)
#need Action's that update the db
# scheduler reads the state, decides which action to do and launches. 
#  the action sends a message to the taskserver to launch a task
# the task runs a script and then updates the db

# I need an Action that does not launch a task, it just sets the current state to Action.state
# taskserver sets the PID
        
class TestScheduler(unittest.TestCase):
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
        s.get_new_active_obs(self.dbi)
        self.assertEqual(len(s.active_obs),self.ntimes*self.npols)
    def test_start(self):
        class SuccessAction(sch.Action):
            def _command(me):
                me.dbi = self.dbi
                me.dbi.set_obs_status(me.task)
    
        def all_done():
            for f in self.files:
                if self.dbi.get_obs_status(f) != 'COMPLETE': return False
        s = sch.Scheduler(nstills=1, actions_per_still=1, blocksize=10)
        t = threading.Thread(target=s.start, args=(self.dbi, SuccessAction))
        t.start()
        tstart = time.time()
        while not all_done() and time.time() - tstart < 1: time.sleep(.1)
        s.quit()
        for f in self.files: self.assertEqual(self.dbi.get_obs_status(f), 'COMPLETE')
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
    def test_clean_completed_actions(self):
        """
        todo
        """
        dbi = FakeDataBaseInterface(10)
        class FakeAction(sch.Action):
            def _command(self):
                dbi.files[self.obs] = self.task
        s = sch.Scheduler(nstills=1, actions_per_still=1, blocksize=10)
        s.get_new_active_obs(self.dbi)
        s.update_action_queue(self.dbi, ActionClass=FakeAction)
        a = s.pop_action_queue(0)
        s.launch_action(a)
        self.assertEqual(len(s.launched_actions[0]), 1)
        s.clean_completed_actions(dbi)
        self.assertEqual(len(s.launched_actions[0]), 0)
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

