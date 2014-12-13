import aipy as a, numpy as n

global data_loc; data_loc = '/Users/mpresley/Research/Research_Adrian_Aaron/gs_data'
global gsm_fits_file; gsm_fits_file = '/Users/mpresley/soft/gsm/haslam408_extrap_fq_0.1_32.fits'
global fig_loc; fig_loc = '/Users/mpresley/soft/capo/mep/figures'

def load_Q_file(gh='grid',del_bl=4.,num_bl=10,beam_sig=0.09,fq=0.05,lmax=10):
    if gh=='grid':
        Q_file = '{0}/Q_matrices/Q_grid_del_bl_{1:.2f}_sqGridSideLen_{2}_beam_sig_{3:.2f}_fq_{4:.3f}.npz'.format(data_loc,del_bl,num_bl,beam_sig,fq)
    elif gh=='hybrid':
        Q_file = '{0}/Q_matrices/Q_hybrid_del_bl_{1:.2f}_sqGridSideLen_{2}_fq_{3:.3f}.npz'.format(data_loc,del_bl,num_bl,fq)
    
    Qstuff = n.load(Q_file)
    Q = Qstuff['Q']
    Q = Q[:,0:(lmax+1)**2]
    lms = Qstuff['lms']
    lms = lms[0:(lmax+1)**2,:]
    baselines = Qstuff['baselines']
    return baselines,Q,lms 

def total_noise_covar(ninst_sig,del_bl=4.,num_bl=10,beam_sig=0.09,fq=0.05):
    Nfg_file = n.load('{0}/gsm_matrices/gsm_grid_del_bl_{1:.2f}_sqGridSideLen_{2}_beam_sig_{3:.2f}_fq_{4:.3f}.npz'.format(data_loc,del_bl,num_bl,beam_sig,fq))
    Nfg = Nfg_file['matrix']
    #Nfg = Nfg[1:,1:]
    #print 'fg ',Nfg.shape
    Ninst = (ninst_sig**2)*n.identity(num_bl*num_bl-1)
    Ntot = Nfg + Ninst 
    return Ntot

def load_mc_data(mc_loc,del_bl=4.,num_bl=10,beam_sig=0.09,fq=0.05):
    """
    the returned ys has dimensions num_bl x num_mc
    """
    mc_mat_loc = '{0}/grid_del_bl_{1:.2f}_num_bl_{2}_beam_sig_{3:.2f}_fq_{4:.3f}'.format(mc_loc,del_bl,num_bl,beam_sig,fq)
    for ii in range(1):
        stuff = n.load('{0}/mc_{1}.npz'.format(mc_mat_loc,ii))
        print 'mc shape ', stuff['matrix'].shape
        if ii==0:
            ys = stuff['matrix']
        else:
            ys = n.hstack((ys,stuff['matrix']))
    return ys 


