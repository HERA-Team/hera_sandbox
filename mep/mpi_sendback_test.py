# MPI wrapper code that acts as a template for one master process
# dishing out tasks to a number of slaves, recieving their finished
# calculations, and giving them another task until all are complete.
#
# The template task is creating entries of a matrix

from mpi4py import MPI 
import numpy as np

# define mpi parameters
comm = MPI.COMM_WORLD
rank = comm.Get_rank() #int associated with each processor. ranges from 0 to number of processors
size=comm.Get_size()
master = 0
num_slaves = size-1
num_sent = 0 
num_messages = num_slaves*2  
message = np.arange(num_messages)*100

# Big running loop
# If I am the master process
if rank==master:
    print "I am the master! Muahaha!"
    # send out first round of messages
    for kk in range(num_slaves):
        print "num_sent = ",num_sent
        comm.send(message[num_sent],dest=kk+1)
        print "message ",message[num_sent]," was sent to slave ",kk+1
        num_sent +=1
    print "Master sent out first round of messages"
    # listen for results and send out new assignments
    for kk in range(num_messages):
        message_back = comm.recv(source=MPI.ANY_SOURCE)
        status = MPI.Status()
        source = status.Get_source()
        print 'Master just received message ',message_back,' from slave ',source
        # if there are more things to do, send out another assignment
        if num_sent<num_messages:
            comm.send(message[num_sent],dest=source)
            print "Master sent out message ",message[num_sent],' to slave ',source
            num_sent +=1
        else:
            # send a -1 to tell slave that task is complete
            comm.send(-1,dest=source)
            print "Master sent out the finished i,j to slave ",source
# If I am a slave and there are not more slaves than jobs
elif rank<=num_messages:
    print "I am slave ",rank
    complete = False
    while not complete:
        # Get message
        recv_message = comm.recv(source=master)
        status = MPI.Status()
        print "slave ",rank," just recieved message ",recv_message
        if recv_message==-1:
            # if there are no more messages
            complete=True
            print "slave ",rank," acknoledges EOM"
        else:
            # send message back
            comm.send(recv_message+1,dest=master)
            print "Slave ",rank," sent back message ",recv_message+1
comm.Barrier()

if rank==master:
    print "The master proclaims completion."

comm.Finalize()

