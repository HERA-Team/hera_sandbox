#! /usr/bin/env python
import ddr_compress as ddr
import logging; logging.basicConfig(level=logging.DEBUG)

#STILLS = ['still0', 'still1', 'still2', 'still3']
#STILLS = ['still3', 'still4', 'still5'] # for RAL test system
#STILLS = ['still4', 'still5','still4','still5'] # for RAL test system
STILLS = ['node18','node22','node22','node22']
PORTS = [14204,14204,14205,14206]
ACTIONS_PER_STILL = 4 # how many actions that run in parallel on a still
BLOCK_SIZE = 10 # number of files that are sent together to a still
TIMEOUT = 600 # seconds; how long a task is allowed to be running before it is assumed to have failed
SLEEPTIME = 1. # seconds; throttle on how often the scheduler polls the database

dbi = ddr.dbi.DataBaseInterface()
task_clients = [ddr.task_server.TaskClient(dbi, s,port=p) for (s,p) in zip(STILLS,PORTS)]
scheduler = ddr.task_server.Scheduler(task_clients, actions_per_still=ACTIONS_PER_STILL, blocksize=BLOCK_SIZE)
scheduler.start(dbi, ActionClass=ddr.task_server.Action, action_args=(task_clients,TIMEOUT), sleeptime=SLEEPTIME)
