print "I exist"
from mpi4py import MPI 
import sys
import aipy as a, numpy as n, pylab as p
import capo as C
import basic_amp_aa_grid_gauss as agg
import useful_functions as uf
import matplotlib as mpl
print "imported everything"

def haslam_extrap(hasdat=None,fq=0.1,save=True):
    alf0=2.8; var = 0.1 # from table on page 4 of http://arxiv.org/pdf/1106.0007.pdf
    if hasdat==None:
        hasmap = a.map.Map(fromfits='/global/homes/m/mpresley/scripts/general_files/fits_files/haslam408_32.fits')
        hasdat = hasmap.map.map 
    alf = n.random.randn(hasdat.shape[0])*var
    print alf 
    fqdat = hasdat*(fq/0.408)**(alf-alf0) 
    hasmap.map.map = fqdat
    if save: hasmap.to_fits('/global/homes/m/mpresley/scripts/general_files/fits_files/haslam408_extrap_fq_{0}_32.fits'.format(fq),clobber=True)
    return fqdat

def generate_sky_model_y(baselines,beamsig,gsm_map=None,gsm_data_file=None):
    """
    y is a vector of the visibilities at different baselines
    """
    if gsm_data_file!=None:
        healmap = a.map.Map(fromfits=gsm_data_file)
    elif gsm_map==None:
        return None
    
    px_array = n.arange(healmap.npix()) # gets an array of healpix pixel indices
    rx,ry,rz = n.array(healmap.px2crd(px_array,ncrd=3)) # finds the topocentric coords for each healpix pixel
    phi,theta = n.array(healmap.px2crd(px_array,ncrd=2)) # phi,theta in math coords
    print px_array.shape
    true_sky = healmap.map.map
    amp = uf.gaussian(beamsig,n.zeros_like(theta),phi)
    dOmega = 4*n.pi/px_array.shape[0]

    visibilities = n.zeros(baselines.shape[0],dtype='complex')
    print baselines.shape[0]
    for kk in range(baselines.shape[0]):
        print kk
        bx,by,bz = baselines[kk]
        Vis = amp*true_sky*n.exp(2j*n.pi*(bx*rx+by*ry+bz*rz))*dOmega
        visibilities[kk] = n.sum(Vis)
    return visibilities


# define mpi parameters
comm = MPI.COMM_WORLD
rank = comm.Get_rank() #int associated with each processor. ranges from 0 to number of processors
size=comm.Get_size()
master = 0
num_slaves = size-1
print "defined mpi params"

# define parameters related to calculation 
hasmap = a.map.Map(fromfits='/global/homes/m/mpresley/scripts/general_files/fits_files/haslam408_32.fits')

_,num0,beam_sig,del_bl,num_bl = sys.argv
beam_sig=float(beam_sig); del_bl=float(del_bl); num_bl=int(num_bl)
baselines = agg.make_pos_array(del_bl,num_bl)

savekey = 'grid_del_bl_{0:.2f}_num_bl_{1}_beam_sig_{2:.2f}_'.format(del_bl,num_bl,beam_sig)

print "defined hasmap stuff"

# define parameters related to task-mastering
num1 = 1
matrix = np.zeros([num0,num1])
assignment_matrix = np.arange(np.prod(matrix.shape)).reshape(matrix.shape)
numToDo = num0*num1
num_sent = 0 # this functions both as a record of how many assignments have 
             # been sent and as a tag marking which matrix entry was calculated

print "defined task-mastering stuff"

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
            gsm_map = haslam_extrap(hasdat=hasmap.map.map,save=False)
            element = generate_sky_model_y(baselines,beamsig,gsm_map=gsm_map)
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



