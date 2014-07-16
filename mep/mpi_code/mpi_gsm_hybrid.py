from mpi4py import MPI 
import sys
import numpy as n
import useful_functions as uf
import aipy as a
import basic_amp_aa_grid_gauss as agg

def compute_element(bli,blj,amp):
    bix,biy,biz = bli; bjx,bjy,bjz = blj
    element = 0
#    print 'amp shape = ',amp.shape
#    print 'px_array shape = ',px_array.shape
    for kk in px_array:
        rx,ry,rz = crd_array[:,kk]          
        Gik = amp[kk]*n.exp(-2j*n.pi*fq*(bix*rx+biy*ry+biz*rz))*dOmega
        Gjk_star = n.conj(amp[kk]*n.exp(-2j*n.pi*fq*(bjx*rx+bjy*ry+bjz*rz)))
        Rkk = Rdata[kk]*Rdata[kk]
        element += Gik*Rkk*Gjk_star
    return element

# define mpi parameters
comm = MPI.COMM_WORLD
rank = comm.Get_rank() #int associated with each processor. ranges from 0 to number of processors
size=comm.Get_size()
master = 0
num_slaves = size-1

# define parameters related to calculation
fq = 0.1
healmap = a.map.Map(fromfits='/global/homes/m/mpresley/scripts/general_files/hi1001_32.fits')
global px_array; px_array = n.arange(healmap.npix()) # gets an array of healpix pixel indices
global crd_array; crd_array = n.array(healmap.px2crd(px_array,ncrd=3)) # finds the topocentric coords for each healpix pixel
global Rdata; Rdata = healmap.map.map
global dOmega; dOmega = 4*n.pi/px_array.shape[0]
phi,theta = n.array(healmap.px2crd(px_array,ncrd=2))
#print 'theta max = ',max(theta)
#print 'phi max = ',max(phi)

maxl = 10
_,del_bl,num_bl = sys.argv
del_bl=float(del_bl);num_bl=int(num_bl)

beamsig_largebm = 10/(2*n.pi*del_bl*(num_bl-1)) #1.0                                         
beamsig_smallbm = 10/(2*n.pi*del_bl) #0.25          
smallbm_inds = (int(n.floor(num_bl/2)),int(n.floor(num_bl/2)))

savekey = 'hybrid_del_bl_{0:.2f}_num_bl_{1}_'.format(del_bl,num_bl)

amp_largebm = uf.gaussian(beamsig_largebm,n.zeros_like(theta),phi)
amp_smallbm = uf.gaussian(beamsig_smallbm,n.zeros_like(theta),phi)
baselines = agg.make_pos_array(del_bl,num_bl)

# define matrix to be calculated
num = len(baselines)
matrix = n.zeros([num,num],dtype=n.complex)
# define parameters related to task-mastering
numToDo = num*(num+1)/2
print 'numToDo = ',numToDo
assn_inds = []
for ii in range(num+1):
    for jj in range(ii+1):
        assn_inds.append((ii,jj))
    num_sent = 0 # this functions both as a record of how many assignments have 
                 # been sent and as a tag marking which matrix entry was calculated
    num_complete = 0
# Big running loop
# If I am the master process
if rank==master:
    print "I am the master! Muahaha!"
    # send out first round of assignments
    for kk in range(num_slaves):
        selectedi, selectedj = assn_inds[kk]
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
        matrix[selectedj,selectedi] = n.conj(entry)
        print 'Master just received element (i,j) = ',selectedi,selectedj,' from slave ',source
        num_complete += 1
        print 'num complete = ',num_complete
        # if there are more things to do, send out another assignment
        if num_sent<numToDo:
            selectedi, selectedj = assn_inds[kk]
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
            if (selectedi,selectedj)==smallbm_inds: amp = amp_smallbm # if it's the baseline with the large smear/small beam                                                               
            else: amp = amp_largebm
            bli = baselines[selectedi,:]
            blj = baselines[selectedj,:]
            element = compute_element(bli,blj,amp)
            # send answer back
            comm.send((rank,element),dest=master)
            comm.send(selectedi,dest=master)
            comm.send(selectedj,dest=master)
            print "Slave ",rank," sent back i,j = ",selectedi,selectedj
comm.Barrier()

if rank==master:
    print "The master will now save the matrix."
    n.savez_compressed('/global/homes/m/mpresley/scripts/gsm_matrices/{0}'.format(savekey),matrix=matrix,baselines=baselines)

MPI.Finalize()

