import SocketServer
import logging, threading, subprocess, time

logger = logging.getLogger('taskserver')

PKT_LINE_LEN = 80

def pad(s, line_len=PKT_LINE_LEN):
    return (s + ' '*line_len)[:line_len]

def to_pkt(task, obs, args):
    nlines = len(args) + 3
    return ''.join(map(pad, [str(nlines), task, obs] + args))

def from_pkt(pkt, line_len=PKT_LINE_LEN):
    nlines,pkt = pkt[:line_len].rstrip(), pkt[line_len:]
    nlines = int(nlines)
    task,pkt = pkt[:line_len].rstrip(), pkt[line_len:]
    obs,pkt = pkt[:line_len].rstrip(), pkt[line_len:]
    args = []
    for i in xrange(nlines-3):
        arg,pkt = pkt[:line_len].rstrip(), pkt[line_len:]
        args.append(arg)
    return task, obs, args

class Task:
    def __init__(self, task, obs, args, dbi, cwd='.'):
        self.task = task
        self.obs = obs
        self.args = args
        self.dbi = dbi
        self.cwd = cwd
        self.process = None
    def run(self):
        if not self.process is None:
            raise RuntimeError('Cannot run a Task that has been run already.')
        self.record_launch()
        self.process = self._run()
    def _run(self):
        return subprocess.Popen(['do_%s.sh' % self.task] + self.args, cwd=self.cwd) # XXX d something with stdout stderr
    def poll(self):
        return self.process.poll()
    def finalize(self):
        self.process.wait()
        if self.poll(): self.record_failure()
        else: self.record_completion()
    def kill(self):
        self.process.kill()
    def record_launch(self): pass # XXX for monitor
    def record_failure(self): pass # XXX for monitor
    def record_completion(self):
        self.dbi.set_obs_status(self.obs, self.task) # XXX change to match actual DBI 
        
'''
    def do_UV(self):
        return self._task('UV',['%s:%s'%(self.pot,self.f)])
    def do_UVC(self):
        return self._task('UVC',[self.f])
    def do_CLEAN_UV(self):
        return self._task('CLEAN_UV',[self.f])
    def do_UVCR(self):
        return self._task('UVCR',[self.f])
    def do_CLEAN_UVC(self):
        return self._task('CLEAN_UVC',[self.f])
    def do_ACQUIRE_NEIGHBORS(self):
        args = []
        if len(self.n1) > 0: args.append(self.n1)
        if len(self.n2) > 0: args.append(self.n2)
        return self._task('ACQUIRE_NEIGHBORS',args) # these will need to have 'still' in filename
    def do_UVCRE(self):
        return self._task('UVCRE',[self.n1, self.f, self.n2])
    def do_NPZ(self):
        return self._task('NPZ',[self.f])
    def do_UVCRR(self):
        return self._task('UVCRR',[self.f])
    def do_NPZ_POT(self):
        return self._task('NPZ_POT',['%s:%s'%(self.pot,self.f)])
    def do_CLEAN_UVCRE(self):
        return self._task('CLEAN_UVCRE',[self.f])
    def do_UVCRRE(self):
        return self._task('UVCRRE',[self.n1, self.f, self.n2])
    def do_CLEAN_UVCRR(self)t
        return self._task('CLEAN_UVCRR',[self.f])
    def do_CLEAN_NPZ(self):
        return self._task('CLEAN_NPZ',[self.f])
    def do_CLEAN_NEIGHBORS(self):
        args = []
        if len(self.n1) > 0: args.append(self.n1)
        if len(self.n2) > 0: args.append(self.n2)
        return self._task('CLEAN_NEIGHBORS',args)
    def do_UVCRRE_POT(self):
        return self._task('UVCRRE_POT',['%s:%s'%(self.pot,self.f)])
    def do_CLEAN_UVCR(self):
        return self._task('CLEAN_UVCR',[self.f])
'''

class TaskHandler(SocketServer.BaseRequestHandler):
    def setup(self):
        logger.info('Connect: %s\n' % str(self.client_address))
    def finish(self):
        logger.info('Disconnect: %s\n' % str(self.client_address))
    def get_pkt(self):
        pkt = self.request.recv(PKT_LINE_LEN)
        nlines = int(pkt)
        pkt += self.request.recv((nlines-1) * PKT_LINE_LEN)
        # XXX should check that pkt is at the correct length
        task, obs, args = from_pkt(pkt)
        return task, obs, args
    def handle(self):
        task, obs, args = self.get_pkt()
        t = Task(task, obs, args, self.server.dbi, self.server.data_dir)
        self.server.append_task(t)
        t.run()
        self.request.send(pad(str(t.pid))) # XXX send back PID for later kill requests?

# XXX swtich to UDP
class TaskServer(SocketServer.TCPServer):
    allow_reuse_address = True
    def __init__(self, dbi, data_dir='.', port=14204, handler=TaskHandler):
        SocketServer.TCPServer.__init__(self, ('', port), handler)
        self.active_tasks_semaphore = threading.Semaphore()
        self.active_tasks = []
        self.dbi = dbi
        self.data_dir = data_dir
        self.is_running = False
    def append_task(self, t):
        self.active_tasks_semaphore.acquire()
        self.active_tasks.append(t)
        self.active_tasks_semaphore.release()
    def finalize_tasks(self, poll_interval=.5):
        while self.is_running:
            self.active_tasks_semaphore.acquire()
            new_active_tasks = []
            for t in self.active_tasks:
                if t.poll() is None: # not complete
                    new_active_tasks.append(t)
                else:
                    t.finalize()
            self.active_tasks = new_active_tasks
            self.active_tasks_semaphore.release()
            time.sleep(poll_interval)
    def kill(self, pid):
        for task in self.active_tasks:
            if task.process.pid == pid:
                task.kill()
                break
    def start(self):
        self.is_running = True
        t = threading.Thread(target=self.finalize_tasks)
        t.start()
        self.serve_forever()
        for task in self.active_tasks: task.kill() # XXX is this cleanup necessary?
        t.join()
    def shutdown(self):
        self.is_running = False
        for t in self.active_tasks:
            try: t.process.kill()
            except(OSError): pass
        SocketServer.TCPServer.shutdown(self)

