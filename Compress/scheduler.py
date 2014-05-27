import time, logging
# XXX deal with 'XXX' in code below
# XXX need to deal with edge files being required on more than one still
# although this isn't necessary for 4 stills, split by polarization

logger = logging.getLogger('scheduler')

FILE_PROCESSING_STAGES  = { # dict linking file status to the next action
    'UV-POT': ('UV-STILL',),
    'UV-STILL': ('UVC',),
    'UVC': ('CLEAN-UV-STILL',),
    'CLEAN-UV-STILL': ('UVCR',),
    'UVCR': ('CLEAN-UVC',),
    'CLEAN-UVC': ('UVCRE',),
    'UVCRE': ('UVCRR',),
    'UVCRR': ('NPZ-POT',),
    'NPZ-POT': ('CLEAN-UVCRE',),
    'CLEAN-UVCRE':, ('UVCRRE-STILL',),
    'UVCRRE-STILL': ('CLEAN-UVCRR',), # do we want uvcRRE to run on UVCR of neighbors, or UVCRR?  I think UVCR, because that's all we get for edge files.
    'CLEAN-UVCRR': ('CLEAN-NPZ-STILL',),
    'CLEAN-NPZ-STILL': ('UVCRRE-POT',),
    'UVCRRE-POT': ('CLEAN-UVCR',),
    'CLEAN-UVCR': ('COMPLETE',),
    'COMPLETE': None,
}

FILE_PROCESSING_PREREQS = { # link task to prerequisite state of neighbors
    'UV-POT': None,
    'UV-STILL': None,
    'UVC': None,
    'CLEAN-UV-STILL': None,
    'UVCR': None,
    'CLEAN-UVC': None,
    'UVCRE': 'UVCR',
    'UVCRR': None,
    'NPZ-POT': None,
    'CLEAN-UVCRE':, None,
    'UVCRRE-STILL': 'UVCR',
    'CLEAN-UVCRR': None,
    'CLEAN-NPZ-STILL': None,
    'UVCRRE-POT': None,
    'CLEAN-UVCR': 'UVCRRE-STILL',
    'COMPLETE': None,
}

def file_to_still(f, nstills):
    '''Return the still that a file should be transferred to.'''
    # XXX
    return

class Action:
    '''An Action performs a task on a file, and is scheduled by a Scheduler.'''
    def __init__(self, f, task, neighbors, still):
        '''f:filename, task:target status, neighbor:adjacent files, 
        still:still action will run on.'''
        self.filename = f
        self.task = task
        self.neighbors = neighbors
        self.still = still
        self.priority = 0
        self.launch_time = -1
    def set_priority(self, p):
        '''Assign a priority to this action.  Highest priorities are scheduled first.'''
        self.priority = p
    def has_prerequisites(self, dbi):
        '''For the given task, check that neighbors are in prerequisite state.
        We don't check that the center file is in the prerequisite state, 
        since this action could not have been generated otherwise.'''
        n_state = FILE_PROCESSING_PREREQS[self.task]
        for n in neighbors:
            status = dbi.get_file_status(n)
            if status != n_state: return False
        return True
    def launch(self):
        '''Run this task.'''
        self.launch_time = time.time()
        # XXX need to actually run this task
        return
    def timed_out(self, timeout=3600.):
        assert(self.launch_time > 0) # Error out if action was not launched
        return time.time() > self.launch_time + timeout
        
        
def action_cmp(x,y): return cmp(x.priority, y.priority)

class Scheduler:
    '''A Scheduler reads a DataBaseInterface to determine what Actions can be
    taken, and then schedules them on stills according to priority.'''
    def __init__(self, nstills=4, actions_per_still=8):
        '''nstills: # of stills in system, 
        actions_per_still: # of actions that can be scheduled simultaneously
                           per still.'''
        self.actions_per_still = actions_per_still
        self.active_files = []
        self._active_file_dict = {}
        self.action_queue = []
        self.launched_actions = {}
        for still in xrange(nstills): self.launched_actions[still] = []
    def start(self, dbi):
        '''Begin scheduling (blocking).
        dbi: DataBaseInterface'''
        logger.info('Beginning scheduler loop')
        while True: # XXX if threading, this should be a self.quit flag
            self.get_new_active_files(dbi)
            self.update_action_queue(dbi)
            # Launch actions that can be scheduled
            for still in self.launched_actions:
                while len(self.launched_actions[still]) < self.actions_per_still:
                    try: a = self.pop_action_queue(still)
                    except(IndexError): # no actions can be taken on this still
                        logger.info('No actions available for still-%d\n' % still)
                        break # move on to next still
                    self.launch_action(a)
            self.clean_completed_actions(dbi)
    def pop_action_queue(self, still):
        '''Return highest priority action for the given still.'''
        return self.launched_actions[still]
    def launch_action(self, a):
        '''Launch the specified Action and record its launch for tracking later.'''
        self.launched_actions[a.still].append(a)
        a.launch()
    def clean_completed_actions(self, dbi):
        '''Check launched actions for completion or timeout.'''
        for still in self.launched_actions:
            updated_actions = []
            for cnt, a in enumerate(self.launched_actions[still]):
                status = dbi.get_file_status(a.filename)
                if status == a.task:
                    logger.info('Task %s for file %s on still %d completed successfully.' % (a.task, a.filename, still)
                    # not adding to updated_actions removes this from list of launched actions
                elif a.timeout(): 
                    logger.info('Task %s for file %s on still %d TIMED OUT.' % (a.task, a.filename, still)
                    # XXX make db entry for documentation
                    # XXX actually kill the process if alive
                else: # still active
                    updated_actions.append(a)
            self.launched_actions[still] = updated_actions
    def already_launched(self, action):
        '''Determine if this action has already been launched.  Enforces
        fact that only one valid action can be taken for a given file
        at any one time.'''
        for a in self.launched_actions[action.still]:
            if a.filename == action.filename: return True
        return False
    def get_new_active_files(self, dbi):
        '''Check for any new files that may have appeared.  Actions for
        these files may potentially take priority over ones currently
        active.'''
        # XXX If actions have been launched since the last time this
        #was called, clean_completed_actions() must be called first to ensure
        #that cleanup occurs before.  Is this true? if so, should add mechanism
        #to ensure ordering
        for f in dbi.ordered_files():
            if not dbi.is_completed(f) and not self._active_file_dict.has_key(f):
                    self._active_file_dict[f] = len(self.active_files)
                    self.active_files.append(f)
    def update_action_queue(self, dbi):
        '''Based on the current list of active files (which you might want
        to update first), generate a prioritized list of actions that 
        can be taken.'''
        actions = [self.get_action(dbi,f) for f in self.active_files]
        actions = [a for a in actions if not a is None] # remove unactionables
        actions = [a for a in actions if self.already_launched(a)] # filter actions already launched
        for i,a in enumerate(actions):
            a.set_priority(self.determine_priority(a,i))
        actions.sort(action_cmp) # place most important actions last
        self.action_queue = actions # completely throw out previous action list
    def get_action(self, dbi, f):
        '''Find the next actionable step for file f (one for which all
        prerequisites have been met.  Return None if no action is available.
        This function is allowed to return actions that have already been
        launched.'''
        status = dbi.get_file_status(f)
        next_step = FILE_PROCESSING_STAGES[status]
        if next_step is None: return None # file is complete
        neighbors = dbi.get_neighbors(f)
        nstills = len(self.launched_actions)
        still = file_to_still(f, nstills) 
        a = Action(f, next_step, neighbors, still)
        if a.has_prerequisites(dbi): return a
        else: return None
    def determine_priority(self, action, fileorder):
        '''Assign a priority to an action based on its status and the time
        order of the file to which this action is attached.'''
        return fileorder # prioritize any possible action on the newest file
        # XXX might want to prioritize finishing a file already started before
        # moving to the latest one (at least, up to a point) to avoid a
        # build up of partial files.  But if you prioritize files already
        # started too excessively, then the queue could eventually fill with
        # partially completed tasks that are failing for some reason
        
                    
                    
                    
