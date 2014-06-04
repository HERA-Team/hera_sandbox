import unittest, subprocess, os, threading, time
import ddr_compress.task_server as ts
import ddr_compress.scheduler as sch

import logging; logging.basicConfig(level=logging.DEBUG)

class FakeDataBaseInterface:
    def __init__(self, nfiles=3):
        self.files = {}
        self.pids = {}
        self.stills = {}
        self.paths = {}
        for i in xrange(nfiles):
            self.files[i] = 'UV_POT'
            self.pids[i] = -1
            self.stills[i] = 'localhost'
            self.paths[i] = os.path.abspath('.')
    def get_obs_status(self, obsnum):
        return self.files[obsnum]
    def list_observations(self):
        files = self.files.keys()
        files.sort()
        return files
    def get_neighbors(self, obsnum):
        n1,n2 = obsnum-1, obsnum+1
        if not self.files.has_key(n1): n1 = None
        if not self.files.has_key(n2): n2 = None
        return (n1,n2)
    def set_obs_status(self, obs, status):
        self.files[obs] = status
    def set_obs_pid(self, obs, pid):
        self.pids[obs] = pid
    def get_obs_pid(self, obs):
        return self.pids[obs]
    def get_input_file(self, obsnum):
        return 'localhost',os.getcwd()+'/data','test%d.uv' % (obsnum+1)
    def get_output_path(self, obsnum):
        return 'localhost','.'
    def get_obs_still_host(self, obsnum):
        return self.stills[obsnum]
    def get_obs_still_path(self, obsnum):
        return self.paths[obsnum]
    def set_obs_still_host(self, obsnum, host):
        self.stills[obsnum] = host
    def set_obs_still_path(self, obsnum, path):
        self.paths[obsnum] = path

class TestTaskScheduler(unittest.TestCase):
    def setUp(self):
        self.dbi = FakeDataBaseInterface()
        self.ts = ts.TaskServer(self.dbi)
        self.ts_thread = threading.Thread(target=self.ts.start)
        self.ts_thread.start()
    def tearDown(self):
        self.ts.shutdown()
        self.ts_thread.join()
        # clean up files made by script chain
        subprocess.call(['rm', '-rf','/tmp/test*.uv*'], shell=True)
        subprocess.call(['rm', '-rf','test*.uv*'], shell=True)
    def test_end_to_end(self):
        tc = {0:ts.TaskClient(self.dbi, 'localhost')}
        s = ts.Scheduler(tc, nstills=1, actions_per_still=1)
        thd = threading.Thread(target=s.start, args=(self.dbi,), kwargs={'ActionClass':ts.Action, 'action_args':(tc,), 'sleeptime':0})
        thd.start()
        try:
            def all_done():
                for f in self.dbi.files:
                    if self.dbi.get_obs_status(f) != 'COMPLETE': return False
                return True
            while not all_done(): time.sleep(.1)
        finally:
            s.quit()
            thd.join()
        for i in self.dbi.files:
            self.assertEqual(self.dbi.get_obs_status(i), 'COMPLETE')
    # XXX need to test killing

if __name__ == '__main__':
    unittest.main()
