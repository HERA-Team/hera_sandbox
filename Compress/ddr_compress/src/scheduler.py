import time, logging

logger = logging.getLogger('scheduler')

FILE_PROCESSING_STAGES = ['NEW','UV_POT', 'UV', 'UVC', 'CLEAN_UV', 'UVCR', 'CLEAN_UVC',
    'ACQUIRE_NEIGHBORS', 'UVCRE', 'NPZ', 'UVCRR', 'NPZ_POT', 'CLEAN_UVCRE', 'UVCRRE',
    'CLEAN_UVCRR', 'CLEAN_NPZ', 'CLEAN_NEIGHBORS', 'UVCRRE_POT', 'CLEAN_UVCR', 'CLEAN_UVCRRE', 'COMPLETE']
FILE_PROCESSING_LINKS = {}
for i,k in enumerate(FILE_PROCESSING_STAGES[:-1]):
    FILE_PROCESSING_LINKS[k] = FILE_PROCESSING_STAGES[i+1]
FILE_PROCESSING_LINKS['COMPLETE'] = None

FILE_PROCESSING_PREREQS = { # link task to prerequisite state of neighbors, key not present assumes no prereqs
    'ACQUIRE-NEIGHBORS': (FILE_PROCESSING_STAGES.index('UVCR'), FILE_PROCESSING_STAGES.index('CLEAN_UVCR')),
    'CLEAN_UVCR': (FILE_PROCESSING_STAGES.index('UVCRRE'),None),
}

class Action:
    '''An Action performs a task on an observation, and is scheduled by a Scheduler.'''
    def __init__(self, f, task, neighbors, still, timeout=3600.):
        '''f:obs, task:target status, neighbor:adjacent obs, 
        still:still action will run on.'''
        self.obs = f
        self.task = task
        self.neighbors = neighbors
        self.still = still
        self.priority = 0
        self.launch_time = -1
        self.timeout = timeout
    def set_priority(self, p):
        '''Assign a priority to this action.  Highest priorities are scheduled first.'''
        self.priority = p
    def has_prerequisites(self, dbi):
        '''For the given task, check that neighbors are in prerequisite state.
        We don't check that the center obs is in the prerequisite state, 
        since this action could not have been generated otherwise.'''
        try: index1,index2 = FILE_PROCESSING_PREREQS[self.task]
        except(KeyError): # this task has no prereqs
            return True
        for n in self.neighbors:
            if n is None: continue # if no neighbor exists, don't wait on it
            status = dbi.get_obs_status(n)
            if status is None: # indicates that obs hasn't been entered into DB yet
                return False
            index = FILE_PROCESSING_STAGES.index(status)
            if not index1 is None and index < index1: return False
            if not index2 is None and index >= index2: return False
        return True
    def launch(self, launch_time=None):
        '''Run this task.'''
        if launch_time is None: launch_time = time.time()
        self.launch_time = launch_time
        return self._command()
    def _command(self):
        '''Replace this function in a subclass to execute different tasks.'''
        return
    def timed_out(self, curtime=None):
        assert(self.launch_time > 0) # Error out if action was not launched
        if curtime is None: curtime = time.time()
        return curtime > self.launch_time + self.timeout
        
def action_cmp(x,y): return cmp(x.priority, y.priority)

class Scheduler:
    '''A Scheduler reads a DataBaseInterface to determine what Actions can be
    taken, and then schedules them on stills according to priority.'''
    def __init__(self, nstills=4, actions_per_still=8, blocksize=10):
        '''nstills: # of stills in system, 
        actions_per_still: # of actions that can be scheduled simultaneously
                           per still.'''
        self.nstills = nstills
        self.actions_per_still = actions_per_still
        self.blocksize = blocksize
        self.active_obs = []
        self._active_obs_dict = {}
        self.action_queue = []
        self.launched_actions = {}
        for still in xrange(nstills): self.launched_actions[still] = []
        self._run = False
    def quit(self):
        self._run = False
    def start(self, dbi, ActionClass=None):
        '''Begin scheduling (blocking).
        dbi: DataBaseInterface'''
        logger.info('Beginning scheduler loop')
        self._run = True
        while self._run:
            tic = time.time()
            self.get_new_active_obs(dbi)
            self.update_action_queue(dbi, ActionClass)
            # Launch actions that can be scheduled
            for still in self.launched_actions:
                while len(self.launched_actions[still]) < self.actions_per_still:
                    try: a = self.pop_action_queue(still)
                    except(IndexError): # no actions can be taken on this still
                        logger.info('No actions available for still-%d\n' % still)
                        break # move on to next still
                    self.launch_action(a)
            self.clean_completed_actions(dbi)
            time.sleep(.1)
    def pop_action_queue(self, still):
        '''Return highest priority action for the given still.'''
        for i in xrange(len(self.action_queue)):
            if self.action_queue[i].still == still:
                return self.action_queue.pop(i)
        raise IndexError('No actions available for still-%d\n' % still)
    def launch_action(self, a):
        '''Launch the specified Action and record its launch for tracking later.'''
        self.launched_actions[a.still].append(a)
        a.launch()
    def clean_completed_actions(self, dbi):
        '''Check launched actions for completion or timeout.'''
        for still in self.launched_actions:
            updated_actions = []
            for cnt, a in enumerate(self.launched_actions[still]):
                status = dbi.get_obs_status(a.obs)
                if status == a.task:
                    logger.info('Task %s for obs %s on still %d completed successfully.' % (a.task, a.obs, still))
                    # not adding to updated_actions removes this from list of launched actions
                elif a.timed_out(): 
                    logger.info('Task %s for obs %s on still %d TIMED OUT.' % (a.task, a.obs, still))
                    # XXX make db entry for documentation
                    # XXX actually kill the process if alive
                else: # still active
                    updated_actions.append(a)
            self.launched_actions[still] = updated_actions
    def already_launched(self, action):
        '''Determine if this action has already been launched.  Enforces
        fact that only one valid action can be taken for a given obs
        at any one time.'''
        for a in self.launched_actions[action.still]:
            if a.obs == action.obs: return True
        return False
    def get_new_active_obs(self, dbi):
        '''Check for any new obs that may have appeared.  Actions for
        these obs may potentially take priority over ones currently
        active.'''
        # XXX If actions have been launched since the last time this
        #was called, clean_completed_actions() must be called first to ensure
        #that cleanup occurs before.  Is this true? if so, should add mechanism
        #to ensure ordering
        for f in dbi.list_observations():
            if dbi.get_obs_status(f) != 'COMPLETE' and not self._active_obs_dict.has_key(f):
                    self._active_obs_dict[f] = len(self.active_obs)
                    self.active_obs.append(f)
    def update_action_queue(self, dbi, ActionClass=None):
        '''Based on the current list of active obs (which you might want
        to update first), generate a prioritized list of actions that 
        can be taken.'''
        actions = [self.get_action(dbi,f,ActionClass=ActionClass) for f in self.active_obs]
        actions = [a for a in actions if not a is None] # remove unactionables
        actions = [a for a in actions if not self.already_launched(a)] # filter actions already launched
        for a in actions: a.set_priority(self.determine_priority(a,dbi))
        actions.sort(action_cmp, reverse=True) # place most important actions first
        self.action_queue = actions # completely throw out previous action list
    def get_action(self, dbi, f, ActionClass=None):
        '''Find the next actionable step for obs f (one for which all
        prerequisites have been met.  Return None if no action is available.
        This function is allowed to return actions that have already been
        launched.
        ActionClass: a subclass of Action, for customizing actions.  
            None defaults to the standard Action'''
        status = dbi.get_obs_status(f) 
        next_step = FILE_PROCESSING_LINKS[status]
        if next_step is None: return None # obs is complete
        neighbors = dbi.get_neighbors(f)
        still = self.obs_to_still(f, dbi) 
        if ActionClass is None: ActionClass = Action
        a = ActionClass(f, next_step, neighbors, still)
        if a.has_prerequisites(dbi): return a
        else: return None
    def determine_priority(self, action, dbi):
        '''Assign a priority to an action based on its status and the time
        order of the obs to which this action is attached.'''
        return action.obs#dbi.get_obs_index(action.obs) # prioritize any possible action on the newest obs
        # XXX might want to prioritize finishing a obs already started before
        # moving to the latest one (at least, up to a point) to avoid a
        # build up of partial obs.  But if you prioritize obs already
        # started too excessively, then the queue could eventually fill with
        # partially completed tasks that are failing for some reason
    def obs_to_still(self, f, dbi):
        '''Return the still that a obs should be transferred to.'''
        cnt = f#dbi.get_obs_index(f)
        return (cnt / self.blocksize) % self.nstills

