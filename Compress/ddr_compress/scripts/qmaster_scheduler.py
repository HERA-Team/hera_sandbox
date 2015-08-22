#! /usr/bin/env python
import ddr_compress as ddr,os,configparser
import logging; logging.basicConfig(level=logging.DEBUG)

#STILLS = ['still0', 'still1', 'still2', 'still3']
#STILLS = ['still3', 'still4', 'still5'] # for RAL test system
#STILLS = ['still4', 'still5','still4','still5'] # for RAL test system
configfile = os.path.expanduser('~/.ddr_compress/still.cfg')
if os.path.exists(configfile):#todo config file an inputable thing above
    config = configparser.ConfigParser()
    configfile = os.path.expanduser(configfile)
    config.read(configfile)
    STILLS = map(str,config['scheduler']['stills'].split(','))
    PORTS = map(int,config['scheduler']['ports'].split(','))
    ACTIONS_PER_STILL = int(config['scheduler']['actions_per_still'])
    BLOCK_SIZE = int(config['scheduler']['block_size'])
    TIMEOUT = int(config['scheduler']['timeout'])
    SLEEPTIME = int(config['scheduler']['sleeptime'])
    print STILLS,PORTS
else:
    STILLS = ['still4', 'still5']
    PORTS = [14204,14204]
    ACTIONS_PER_STILL = 8 # how many actions that run in parallel on a still
    BLOCK_SIZE = 10 # number of files that are sent together to a still
    TIMEOUT = 600 # seconds; how long a task is allowed to be running before it is assumed to have failed
    SLEEPTIME = 1. # seconds; throttle on how often the scheduler polls the database

dbi = ddr.dbi.DataBaseInterface()
task_clients = [ddr.task_server.TaskClient(dbi, s,port=p) for (s,p) in zip(STILLS,PORTS)]
scheduler = ddr.task_server.Scheduler(task_clients, actions_per_still=ACTIONS_PER_STILL,blocksize=BLOCK_SIZE,nstills=len(STILLS))
scheduler.start(dbi, ActionClass=ddr.task_server.Action, action_args=(task_clients,TIMEOUT), sleeptime=SLEEPTIME)
