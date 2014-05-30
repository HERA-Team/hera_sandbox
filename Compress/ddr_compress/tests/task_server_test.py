import unittest, threading
import subprocess, os, time, socket
import ddr_compress.task_server as ts

class SleepTask(ts.Task):
    def _run(self):
        return subprocess.Popen(['sleep','100'], stdout=open(os.devnull,'w'))

class NullTask(ts.Task):
    def _run(self):
        return subprocess.Popen(['ls'], stdout=open(os.devnull,'w'), cwd=self.cwd)

class FakeDataBaseInterface:
    def __init__(self, nfiles=10):
        self.files = {}
        for i in xrange(nfiles):
            self.files[str(i)] = 'UV-POT'
    def get_obs_status(self, filename):
        return self.files[filename]
    def get_obs_index(self, filename):
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
    def set_obs_status(self, filename, status):
        self.files[filename] = status

class TestFunctions(unittest.TestCase):
    def test_pad(self):
        self.assertEqual(len(ts.pad('', 80)), 80)
        self.assertEqual(len(ts.pad('abc'*10, 30)), 30)
    def test_to_pkt(self):
        pkt = ts.to_pkt('UV', 'filename', ['1','2','3'])
        self.assertEqual(len(pkt), 6*ts.PKT_LINE_LEN)
        self.assertEqual(pkt[:ts.PKT_LINE_LEN], ts.pad('6'))
    def test_from_pkt(self):
        pkt = ts.pad('4') + ts.pad('UV') + ts.pad('filename') + ts.pad('1')
        task, obs, args = ts.from_pkt(pkt)
        self.assertEqual(task, 'UV')
        self.assertEqual(obs, 'filename')
        self.assertEqual(args, ['1'])
    def test_to_from_pkt(self):
        pkt = ts.to_pkt('UV', 'filename', ['1','2','3'])
        task, obs, args = ts.from_pkt(pkt)
        self.assertEqual(task, 'UV')
        self.assertEqual(obs, 'filename')
        self.assertEqual(args, ['1','2','3'])
        
class TestTask(unittest.TestCase):
    def setUp(self):
        self.var = 0
        class VarTask(ts.Task):
            def _run(me):
                self.var += 1
                return subprocess.Popen(['ls'], stdout=open(os.devnull,'w'), 
                    cwd=me.cwd)
        self.VarTask = VarTask
    def test_run(self):
        dbi = FakeDataBaseInterface()
        t = self.VarTask('UV', '1', ['filename'], dbi)
        self.assertEqual(t.process, None)
        var = self.var
        t.run()
        self.assertEqual(self.var, var+1)
        self.assertTrue(type(t.process) is subprocess.Popen)
        t.finalize()
        self.assertEqual(dbi.get_obs_status('1'), 'UV')
        self.assertRaises(RuntimeError, t.run)
    def test_kill(self):
        dbi = FakeDataBaseInterface()
        t = SleepTask('UV','1',[],dbi)
        start_t = time.time()
        t.run()
        t.kill()
        t.finalize()
        end_t = time.time()
        self.assertEqual(t.poll(), -9)
        self.assertLess(end_t-start_t, 100)
        self.assertEqual(dbi.get_obs_status('1'), 'UV-POT')
        
class TestTaskServer(unittest.TestCase):
    def setUp(self):
        self.dbi = FakeDataBaseInterface()
    def test_basics(self):
        s = ts.TaskServer(self.dbi)
        t = SleepTask('UV','1',[],self.dbi)
        s.append_task(t)
        self.assertEqual(len(s.active_tasks), 1)
        t.run()
        s.kill(t.process.pid)
        while t.poll() is None: time.sleep(.01)
        self.assertEqual(t.poll(), -9)
        thd = threading.Thread(target=s.finalize_tasks, args=(.1,))
        s.is_running = True
        thd.start()
        s.is_running = False
        thd.join()
        self.assertEqual(len(s.active_tasks),0)
    def test_shutdown(self):
        s = ts.TaskServer(self.dbi)
        t = threading.Thread(target=s.start)
        t.start()
        s.shutdown()
        t.join()
        self.assertFalse(s.is_running)
    def test_send_task(self):
        self.var = 0
        class SleepHandler(ts.TaskHandler):
            def handle(me):
                self.var += 1
                t = SleepTask('UV','1',[],self.dbi)
                t.run()
                me.server.append_task(t)
        s = ts.TaskServer(self.dbi, handler=SleepHandler)
        thd = threading.Thread(target=s.start)
        thd.start()
        try:
            self.assertEqual(len(s.active_tasks), 0)
            self.assertEqual(self.var, 0)
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.connect(('localhost', 14204))
            sock.send('test')
            sock.close()
            while self.var != 1: time.sleep(.1)
            self.assertEqual(self.var, 1)
            self.assertEqual(len(s.active_tasks), 1)
        finally:
            s.shutdown()
            thd.join()
    def test_dbi(self):
        self.var = 0
        for f in self.dbi.files: self.dbi.files[f] = 'UV-POT'
        class NullHandler(ts.TaskHandler):
            def handle(me):
                task, obs, args = me.get_pkt()
                t = NullTask(task, obs, args, self.dbi)
                me.server.append_task(t)
                t.run()
                self.var += 1
        s = ts.TaskServer(self.dbi, handler=NullHandler)
        thd = threading.Thread(target=s.start)
        thd.start()
        try:
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.connect(('localhost', 14204))
            sock.send(ts.to_pkt('UV','1',[]))
            sock.close()
            while self.var != 1: time.sleep(.6)
            self.assertEqual(self.var, 1)
            self.assertEqual(self.dbi.get_obs_status('1'), 'UV')
        finally:
            s.shutdown()
            thd.join()
                

        

if __name__ == '__main__':
    unittest.main()
