# MPI wrapper code that acts as a template for one master process
# dishing out tasks to a number of slaves, recieving their finished
# calculations, and giving them another task until all are complete.
#
# The template task is creating entries of a matrix

from mpi4py import MPI 
import numpy as np


def compute_element(ii,jj):
    return ii*ii+jj*jj

# define mpi parameters
comm = MPI.COMM_WORLD
rank = comm.Get_rank() #int associated with each processor. ranges from 0 to number of processors
size=comm.Get_size()
master = 0
num_slaves = size-1

# define parameters related to calculation 
num0,num1 = 5,5
matrix = np.zeros([num0,num1])
assignment_matrix = np.arange(np.prod(matrix.shape)).reshape(matrix.shape)

# define parameters related to task-mastering
numToDo = num0*num1
num_sent = 0 # this functions both as a record of how many assignments have 
             # been sent and as a tag marking which matrix entry was calculated
    

# Big running loop
# If I am the master process
if rank==master:
    print "I am the master! Muahaha!"
    # send out first round of assignments
    for kk in range(num_slaves):
        selectedi, selectedj = np.where(assignment_matrix==num_sent)
        selectedi = selectedi[0]; selectedj = selectedj[0]
        print "num_sent = ",num_sent
        comm.send(selectedi,dest=kk+1)
        comm.send(selectedj,dest=kk+1)
        print "i,j = ",selectedi,selectedj," was sent to slave ",kk+1
        num_sent +=1
    print "Master sent out first round of assignments"
    # listen for results and send out new assignments
    for kk in range(numToDo):
        source,entry = comm.recv(source=MPI.ANY_SOURCE)
        selectedi = comm.recv(source=source)
        selectedj = comm.recv(source=source)
        # stick entry into matrix 
        matrix[selectedi,selectedj] = entry
        print 'Master just received element (i,j) = ',selectedi,selectedj,' from slave ',source
        # if there are more things to do, send out another assignment
        if num_sent<numToDo:
            selectedi, selectedj = np.where(assignment_matrix==num_sent)
            selectedi = selectedi[0]; selectedj = selectedj[0]
            comm.send(selectedi,dest=source)
            comm.send(selectedj,dest=source)
            print "Master sent out i,j = ",selectedi, selectedj,' to slave ',source
            num_sent +=1
        else:
            # send a -1 to tell slave that task is complete
            comm.send(-1,dest=source)
            comm.send(-1,dest=source)
            print "Master sent out the finished i,j to slave ",source
# If I am a slave and there are not more slaves than jobs
elif rank<=numToDo:
    print "I am slave ",rank
    complete = False
    while not complete:
        # Get assignment
        selectedi = comm.recv(source=master)
        selectedj = comm.recv(source=master)
        print "slave ",rank," just recieved i,j = ",selectedi,selectedj
        if selectedi==-1:
            # if there are no more jobs
            complete=True
            print "slave ",rank," acknoledges job completion"
        else:
            # compute the matrix element
            element = compute_element(selectedi,selectedj)
            # send answer back
            comm.send((rank,element),dest=master)
            comm.send(selectedi,dest=master)
            comm.send(selectedj,dest=master)
            print "Slave ",rank," sent back i,j = ",selectedi,selectedj
comm.Barrier()

if rank==master:
    print "The master will now save the matrix."
    np.savez_compressed('./test_matrix',matrix=matrix)
    print matrix

MPI.Finalize()

