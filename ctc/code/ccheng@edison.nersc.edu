#!/usr/bin/env python

from mpi4py import MPI
import numpy

#set up MPI

comm = MPI.COMM_WORLD #get MPI communicator object
size = comm.size      #total number of processors
rank = comm.rank      #rank of a process
status = MPI.Status() #MPI status object (contains source and tag)

#start of code

if rank == 0:
  
    #integration parameters
    lower_bound = 0
    upper_bound = 5
    step_size = 1
    int_range = numpy.arange(lower_bound,upper_bound,step_size)

    tasks = int_range
    final = numpy.zeros_like(tasks)
    task_index = 0
    num_workers = size-1
    closed_workers = 0

    print "Master starting with %d workers and %d tasks" % (num_workers, len(tasks))

    while closed_workers < num_workers:
        worker_data = comm.recv(source=MPI.ANY_SOURCE,tag=MPI.ANY_TAG,status=status)
        tag = status.Get_tag()
        source = status.Get_source()
        if tag == 0: #worker is ready
            if task_index < len(tasks):
                print "Sending task %d to worker %d" % (task_index,source)
                comm.send((tasks[task_index],(task_index)),dest=source,tag=3)
                task_index += 1
            else:
                comm.send(None,dest=source,tag=2)
        elif tag == 1: #done tag
            results = worker_data
            final[results[1]] = results[0]
            print "Got data from worker %d" % source
        elif tag == 2: #no more workers needed tag
            print "Worker %d exited" % source
            closed_workers +=1 
            
    print "Master finishing" #final computation
    print "Answer = ",numpy.sum(final)*step_size

else:
    
    print "I am a worker with rank %d" % rank
    while True:      
        comm.send(None,dest=0,tag=0) #worker stops until message is received
        task = comm.recv(source=0,tag=MPI.ANY_TAG,status=status)
        tag = status.Get_tag()
        if tag == 3: #start tag
            #do work here
            result = task[0]**2
            comm.send((result,task[1]),dest=0,tag=1) #done tag
        elif tag == 2:
            break
        
    comm.send(None,dest=0,tag=2)
    
