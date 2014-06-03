import SocketServer
import logging, threading, subprocess, time
import socket

logger = logging.getLogger('taskserver')

PKT_LINE_LEN = 80
STILL_PORT = 14204

def pad(s, line_len=PKT_LINE_LEN):
    return (s + ' '*line_len)[:line_len]

def to_pkt(task, obs, args):
    nlines = len(args) + 3
    return ''.join(map(pad, [str(nlines), task, str(obs)] + args))

def from_pkt(pkt, line_len=PKT_LINE_LEN):
    nlines,pkt = pkt[:line_len].rstrip(), pkt[line_len:]
    nlines = int(nlines)
    task,pkt = pkt[:line_len].rstrip(), pkt[line_len:]
    obs,pkt = int(pkt[:line_len].rstrip()), pkt[line_len:]
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
        self.process = self._run()
        self.record_launch()
    def _run(self):
        logger.info('Task._run: (%s,%d) %s' % (self.task,self.obs,' '.join(['do_%s.sh' % self.task] + self.args)))
        return subprocess.Popen(['do_%s.sh' % self.task] + self.args, cwd=self.cwd) # XXX d something with stdout stderr
    def poll(self):
        return self.process.poll()
    def finalize(self):
        self.process.wait()
        if self.poll(): self.record_failure()
        else: self.record_completion()
    def kill(self):
        self.record_failure()
        self.process.kill()
    def record_launch(self): 
        self.dbi.set_obs_pid(self.obs, self.process.pid)
    def record_failure(self):
        self.dbi.set_obs_pid(self.obs, -1)
    def record_completion(self):
        self.dbi.set_obs_status(self.obs, self.task)

class TaskClient:
    def __init__(self, dbi, host, port=STILL_PORT):
        self.dbi = dbi
        self.host_port = (host,port)
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    def _tx(self, task, obs, args):
        pkt = to_pkt(task, obs, args)
        self.sock.sendto(pkt, self.host_port)
    def gen_args(self, task, obs):
        pot,path,basename = self.dbi.get_input_file(obs)
        outhost,outpath = self.dbi.get_output_path(obs)
        stillhost,stillpath = self.dbi.get_still_host(obs), self.dbi.get_still_path(obs)
        neighbors = [(self.dbi.get_still_host(n),self.dbi.get_still_path(n)) + self.dbi.get_input_file(n)
            for n in self.dbi.get_neighbors(obs) if not n is None]
        neighbors_base = list(self.dbi.get_neighbors(obs))
        if not neighbors_base[0] is None: neighbors_base[0] = self.dbi.get_input_file(neighbors_base[0])[-1]
        if not neighbors_base[1] is None: neighbors_base[1] = self.dbi.get_input_file(neighbors_base[1])[-1]
        def interleave(filename, appendage='cR'):
            rv = [filename]
            if not neighbors_base[0] is None: rv = [neighbors_base[0]+appendage] + rv
            if not neighbors_base[1] is None: rv = rv + [neighbors_base[1]+appendage]
            return rv
        args = {
            'UV': [basename, '%s:%s/%s' % (pot,path,basename)],
            'UVC': [basename],
            'CLEAN_UV': [basename],
            'UVCR': [basename+'c'],
            'CLEAN_UVC': [basename+'c'],
            'ACQUIRE_NEIGHBORS': ['%s:%s/%s' % (n[0], n[1], n[-1]+'cR') for n in neighbors if n[0] != stillhost],
            'UVCRE': interleave(basename+'cR'),
            'NPZ': [basename+'cRE'],
            'UVCRR': [basename+'cR'],
            'NPZ_POT': [basename+'cRE.npz', '%s:%s' % (outhost,outpath)],
            'CLEAN_UVCRE': [basename+'cRE'],
            'UVCRRE': interleave(basename+'cRR'),
            'CLEAN_UVCRR': [basename+'cRR'],
            'CLEAN_NPZ': [basename+'cRE.npz'],
            'CLEAN_NEIGHBORS': [n[-1]+'cR' for n in neighbors if n[0] != stillhost],
            'UVCRRE_POT': [basename+'cRRE', '%s:%s' % (outhost,outpath)],
            'CLEAN_UVCR': [basename+'cR'],
            'CLEAN_UVCRRE': [basename+'cRRE'],
            'COMPLETE': [],
        }
        return args[task]
    def tx(self, task, obs):
        args = self.gen_args(task, obs)
        self._tx(task, obs, args)
    def tx_kill(self, obs):
        self._tx('KILL', obs, [self.dbi.get_obs_pid(obs)])

# XXX consider moving this class to a separate file
import scheduler
class Action(scheduler.Action):
    def __init__(self, obs, task, neighbor_status, still, task_client, timeout=3600.):
        scheduler.Action.__init__(self, obs, task, neighbor_status, still, timeout=timeout)
        self.task_client = task_client
    def _command(self):
        logger.debug('Action: task_client(%s,%d)' % (self.task, self.obs))
        self.task_client.tx(self.task, self.obs)

class TaskHandler(SocketServer.BaseRequestHandler):
    def setup(self):
        #logger.debug('Connect: %s\n' % str(self.client_address))
        return
    def finish(self):
        #logger.debug('Disconnect: %s\n' % str(self.client_address))
        return
    def get_pkt(self):
        pkt = self.request[0]
        task, obs, args = from_pkt(pkt)
        return task, obs, args
    def handle(self):
        task, obs, args = self.get_pkt()
        if task == 'KILL':
            self.server.kill(args[0])
            return
        elif task == 'COMPLETE':
            self.server.dbi.set_obs_status(obs, task)
            return
        t = Task(task, obs, args, self.server.dbi, self.server.data_dir)
        self.server.append_task(t)
        t.run()

class TaskServer(SocketServer.UDPServer):
    allow_reuse_address = True
    def __init__(self, dbi, data_dir='.', port=STILL_PORT, handler=TaskHandler):
        SocketServer.UDPServer.__init__(self, ('', port), handler)
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
        SocketServer.UDPServer.shutdown(self)

