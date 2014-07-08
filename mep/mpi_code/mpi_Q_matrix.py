from mpi4py import MPI 
import aipy as a, numpy as n
import capo as C
import useful_functions as uf
from scipy import special

def get_baselines(aa, nb=None):
    """
    This function takes in an antenna array, and the number of baselines 
    to return (returns all if nb=None).
    """
    na = len(aa.ants) # number of antennas
    baselines = n.zeros([(na*na-na)/2,4])
    ll=0
    for ii in n.arange(na):
        for jj in n.arange(ii):
            bx,by,bz = aa.get_baseline(ii,jj,'r') 
            #the baseline for antennas ii and jj 
            bb = n.sqrt(bx*bx+by*by+bz*bz)
            #print ii,jj,[bx,by,bz,bb]
            baselines[ll] = [bx,by,bz,bb]
            ll+=1
    sorted_baselines = baselines[baselines[:,-1].argsort()] # sorts by last column
    sorted_baselines = sorted_baselines[:,0:3]
    if nb!=None: sorted_baselines = sorted_baselines[0:nb,:]
    return sorted_baselines

def get_single_Q_element(im,(tx,ty,tz),amp,baseline,l,m):
    bx,by,bz = baseline
    # compute spherical harmonic
    Ynorm = special.sph_harm(0,0,0,0)
    tx,ty,tz = tx.flatten(),ty.flatten(),tz.flatten()
    theta = n.arctan(ty/tx) # using math convention of theta=[0,2pi], phi=[0,pi]
    phi = n.arccos(n.sqrt(1-tx*tx-ty*ty))
    Y = n.array(special.sph_harm(m,l,theta,phi))/Ynorm #using math convention of theta=[0,2pi], phi=[0,pi]    
    #fringe pattern
    phs = n.exp(-2j*n.pi*fq * (bx*tx+by*ty+bz*tz)) 
    # not sure if I actually need this stuff
    tx.shape = im.uv.shape
    valid = n.logical_not(tx.mask)
    tx = tx.flatten()
    Y.shape = phs.shape = amp.shape = im.uv.shape
    amp = n.where(valid, amp, n.zeros_like(amp))
    phs = n.where(valid, phs, n.zeros_like(phs))
    Y = n.where(valid, Y, n.zeros_like(Y)) 
    Q_element = n.sum(amp*Y*phs)/n.sum(amp)
    return Q_element

# define mpi parameters
comm = MPI.COMM_WORLD
rank = comm.Get_rank() #int associated with each processor. ranges from 0 to number of processors
size=comm.Get_size()
master = 0
num_slaves = size-1

# define parameters related to calculation 
maxl = 15
calfile = 'basic_amp_aa_grid_gauss'
beamsig = 10.*n.pi/180.
savekey = 'grid_5_ns_gauss_10_deg'
fq = 0.1

aa = a.cal.get_aa(calfile, n.array([.10])) #get antenna array
im = a.img.Img(size=200, res=.5) #make an image of the sky to get sky coords
tx,ty,tz = im.get_top(center=(200,200)) #get coords of the zenith?
valid = n.logical_not(tx.mask)
tx,ty,tz = tx.flatten(),ty.flatten(),tz.flatten()
theta = n.arctan(ty/tx) # using math convention of theta=[0,2pi], phi=[0,pi]
phi = n.arccos(n.sqrt(1-tx*tx-ty*ty))
amp = uf.gaussian(beamsig,n.zeros_like(theta),phi) 

baselines = get_baselines(aa)
num0,num1 = len(baselines),(maxl+1)*(maxl+1)
print "num baselines = {0}\n num lms = {1}".format(num0,num1)
lms = n.zeros([num1,2])
ii=0
for ll in range(maxl+1):
	for mm in range(-ll,ll+1):
		lms[ii] = n.array([ll,mm])
		ii+=1
matrix = n.zeros([num0,num1],dtype=n.complex)
assignment_matrix = n.arange(n.prod(matrix.shape)).reshape(matrix.shape)

# compute normalization for spherical harmonics
Ynorm = special.sph_harm(0,0,0,0)

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
        selectedi, selectedj = n.where(assignment_matrix==num_sent)
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
            selectedi, selectedj = n.where(assignment_matrix==num_sent)
            selectedi = selectedi[0]; selectedj = selectedj[0]
            comm.send(selectedi,dest=source)
            comm.send(selectedj,dest=source)
            #print "Master sent out i,j = ",selectedi, selectedj,' to slave ',source
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
            element = get_single_Q_element(im,(tx,ty,tz),amp,baselines[selectedi],lms[selectedj,0],lms[selectedj,1])
            # send answer back
            comm.send((rank,element),dest=master)
            comm.send(selectedi,dest=master)
            comm.send(selectedj,dest=master)
            print "Slave ",rank," sent back i,j = ",selectedi,selectedj
comm.Barrier()

if rank==master:
    n.savez_compressed('./Q_matrices/{0}Q_max_l_{1}'.format(savekey,maxl),Q=matrix,baselines=baselines,lms=lms)
    print "The master has saved the matrix."

MPI.Finalize()


