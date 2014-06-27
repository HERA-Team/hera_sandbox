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

# Big running loop
# If I am the master process
if rank==master:
    num0,num1 = 5,5
    # define parameters related to calculation 
    matrix = np.zeros([num0,num1])
    assignment_matrix = np.arange(np.prod(matrix.shape)).reshape(matrix.shape)
    
    # define parameters related to task-mastering
    numToDo = num0*num1
    num_sent = 0 # this functions both as a record of how many assignments have 
                 # been sent and as a tag marking which matrix entry was calculated

    # send out first round of assignments
    for kk in range(num_slaves):
        selectedi, selectedj = np.where(assignment_matrix==num_sent)
        selectedi = selectedi[0]; selectedj = selectedj[0]
        comm.send(selectedi,dest=kk+1,tag=num_sent)
        comm.send(selectedj,dest=kk+1,tag=num_sent)
        num_sent +=1

    # listen for results and send out new assignments
    for kk in range(numToDo):
        entry = comm.recv(source=MPI.ANY_SOURCE)
        source = MPI.Status.Get_source()
        tag = MPI.Status.Get_tag()
        # figure out entry i,j based on tag
        selectedi, selectedj = np.where(assignment_matrix==tag)
        selectedi = selectedi[0]; selectedj = selectedj[0]
        # stick entry into matrix 
        matrix[selectedi,selectedj] = entry
        print 'Just finished element (i,j) = ',selectedi,selectedj

        # if there are more things to do, send out another assignment
        if num_sent<numToDo:
            selectedi, selectedj = np.where(assignment_matrix==num_sent)
            selectedi = selectedi[0]; selectedj = selectedj[0]
            comm.send(selectedi,dest=source,tag=num_sent)
            comm.send(selectedj,dest=source,tag=num_sent)
            num_sent +=1
        else:
            # send a -1 to tell slave that task is complete
            comm.send(-1,dest=source)
            comm.send(-1,dest=source)
# If I am a slave and there are not more slaves than jobs
elif rank<=numToDo:
    complete = False
    while not complete:
        # Get assignment
        comm.recv(selectedi,source=master,tag=MPI.ANY_TAG)
        comm.recv(selectedj,source=master,tag=MPI.ANY_TAG)
        tag = MPI.Status.Get_tag()
        if selectedi==-1:
            # if there are no more jobs
            complete=True
        else:
            # compute the matrix element
            element = compute_element(selectedi,selectedj)
            # send answer back
            comm.send(element,dest=master,tag=tag)


comm.Barrier()

if rank==master:
    np.savez_compressed(savepath,matrix=matrix)

comm.Finalize()


def compute_element(ii,jj):
    return ii*ii+jj*jj
