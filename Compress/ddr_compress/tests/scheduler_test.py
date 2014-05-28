import unittest
import ddr_compress.scheduler as sch
import random

class NullAction(sch.Action):
    def _command(self): return

class FakeDataBaseInterface:
    def __init__(self, nfiles=10):
        self.files = {}
        for i in xrange(nfiles):
            self.files[str(i)] = 'UV-POT'
    def get_file_status(self, filename):
        return self.files[filename]
    def file_index(self, filename):
        return int(filename)
    def ordered_files(self):
        files = self.files.keys()
        files.sort()
        return files
    def is_completed(self, filename):
        return self.files[filename] == 'COMPLETE'
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
        a = NullAction(self.files[1], self.task, [self.files[0],self.files[2]], self.still)
        self.assertRaises(AssertionError, a.timed_out)
        t0 = 1000
        a.launch(launch_time=t0)
        self.assertFalse(a.timed_out(timeout=100, curtime=t0))
        self.assertTrue(a.timed_out(timeout=100, curtime=t0+110))
    def test_action_cmp(self):
        priorities = range(10)
        actions = [sch.Action(self.files[1], self.task, [self.files[0],self.files[2]], self.still) for p in priorities]
        random.shuffle(priorities)
        for a,p in zip(actions,priorities): a.set_priority(p)
        actions.sort(cmp=sch.action_cmp)
        for cnt,a in enumerate(actions):
            self.assertEqual(a.priority, cnt)
        
class TestScheduler(unittest.TestCase):
    def setUp(self):
        self.nfiles = 10
        dbi = FakeDataBaseInterface(self.nfiles)
        class FakeAction(sch.Action):
            def _command(self):
                dbi.files[self.filename] = self.task
        self.FakeAction = FakeAction
        self.dbi = dbi
    def test_attributes(self):
        s = sch.Scheduler(nstills=1, actions_per_still=1)
        self.assertEqual(s.launched_actions.keys(), [0])
    def test_get_new_active_files(self):
        s = sch.Scheduler(nstills=1, actions_per_still=1)
        s.get_new_active_files(self.dbi)
        for i in xrange(self.nfiles):
            self.assertTrue(str(i) in s.active_files)
    def test_get_action(self):
        s = sch.Scheduler(nstills=1, actions_per_still=1)
        f = '1'
        a = s.get_action(self.dbi, f, ActionClass=self.FakeAction)
        self.assertNotEqual(a, None) # everything is actionable in this test
        self.assertEqual(a.task, sch.FILE_PROCESSING_STAGES[self.dbi.files[f]][0]) # check this links to the next step
    def test_update_action_queue(self):
        s = sch.Scheduler(nstills=1, actions_per_still=1, blocksize=10)
        s.get_new_active_files(self.dbi)
        s.update_action_queue(self.dbi)
        self.assertEqual(len(s.action_queue), self.nfiles)
        self.assertGreater(s.action_queue[0].priority, s.action_queue[-1].priority)
        for a in s.action_queue: self.assertEqual(a.task, 'UV')

if __name__ == '__main__':
    unittest.main()
