import aipy as a, numpy as n, pylab as p
import capo as C
import useful_functions as uf
import global_sky_model as gsm
import sph_harm_coeffs as shc 
import matplotlib as mpl
import healpy as hp

def total_noise_covar(ninst_sig,num_bl,fg_file):
    Nfg_file = n.load(fg_file)
    Nfg = Nfg_file['matrix']
    #Nfg = Nfg[1:,1:]
    #print 'fg ',Nfg.shape
    Ninst = (ninst_sig**2)*n.identity(num_bl)
    Ntot = Nfg + Ninst 
    return Ntot

def error_covariance(Q,N):
    Q = n.matrix(Q)
    N = n.matrix(N)
    Ninv = uf.pseudo_inverse(N,num_remov=N.shape[0]-5)
    info = n.dot(Q.H,n.dot(Ninv,Q))
    err_cov = uf.pseudo_inverse(info)
    err_cov = n.array(err_cov)
    return err_cov

def Q_N_plots(Q,N,lms,save_tag):
    ls = n.arange(max(lms[:,0])+1)
    l_locs = ls*ls

    Q = n.matrix(Q)
    N = n.matrix(N)
    Ninv = uf.pseudo_inverse(N,num_remov=N.shape[0]-5)
    foo = n.dot(Q.H,Q)
    foo = n.array(uf.pseudo_inverse(foo))
    bar = n.dot(Q.H,n.dot(Ninv,Q))
    bar_inv = uf.pseudo_inverse(bar)
    sigma = n.dot(bar_inv,n.dot(bar,bar_inv.H))
    bar_inv = n.array(bar_inv); sigma = n.array(sigma)

    p.imshow(n.log(n.absolute(foo)))
    p.title('(Qdag Q)^-1')
    p.xticks(l_locs,ls)
    p.yticks(l_locs,ls)
    p.xlabel('l')
    p.ylabel('l')
    p.savefig('./figures/{0}_QdagQ_inv.pdf'.format(save_tag))
    p.clf()

    p.scatter(lms[:,0],n.diag(foo),c=lms[:,1],cmap=mpl.cm.PiYG,s=50)
    p.xlabel('l (color is m)')
    p.ylabel('diagonal elements of (Qdag Q)^-1')
    p.colorbar()
    p.savefig('./figures/{0}_QdagQ_inv_diagonal.pdf'.format(save_tag))
    p.clf()

    p.imshow(n.log(n.absolute(bar_inv)))
    p.title('(Qdag Ninv Q)^-1')
    p.xticks(l_locs,ls)
    p.yticks(l_locs,ls)
    p.xlabel('l')
    p.ylabel('l')
    p.savefig('./figures/{0}_QdagNinvQ_inv.pdf'.format(save_tag))
    p.clf()

    p.scatter(lms[:,0],n.diag(bar_inv),c=lms[:,1],cmap=mpl.cm.PiYG,s=50)
    p.xlabel('l (color is m)')
    p.ylabel('diagonal elements of (Qdag Ninv Q)^-1')
    p.colorbar()
    p.savefig('./figures/{0}_QdagNinvQ_inv_diagonal.pdf'.format(save_tag))
    p.clf()

    p.imshow(n.log(n.absolute(sigma)))
    p.title('Sigma Matrix')
    p.xticks(l_locs,ls)
    p.yticks(l_locs,ls)
    p.xlabel('l')
    p.ylabel('l')
    p.savefig('./figures/{0}_sigma.pdf'.format(save_tag))
    p.clf()

    p.scatter(lms[:,0],n.diag(sigma),c=lms[:,1],cmap=mpl.cm.PiYG,s=50)
    p.xlabel('l (color is m)')
    p.ylabel('diagonal elements of Sigma Matrix')
    p.colorbar()
    p.savefig('./figures/{0}_sigma_diagonal.pdf'.format(save_tag))
    p.clf()

def haslam_extrap(hasdat=None,fq=0.1,save=True):
    alf0=2.8; var = 0.1 # from table on page 4 of http://arxiv.org/pdf/1106.0007.pdf
    if hasdat==None:
        hasmap = a.map.Map(fromfits='/Users/mpresley/soft/gsm/haslam408_32.fits')
        hasdat = hasmap.map.map 
    alf = n.random.randn(hasdat.shape[0])*var
    #print alf 
    fqdat = hasdat*(fq/0.408)**(alf-alf0) 
    hasmap.map.map = fqdat
    if save: hasmap.to_fits('/Users/mpresley/soft/gsm/haslam408_extrap_fq_{0}_32.fits'.format(fq),clobber=True)
    return fqdat

def generate_sky_model_y(baselines,beamsig,gsm_data_file):
    """
    y is a vector of the visibilities at different baselines
    """
    healmap = a.map.Map(fromfits=gsm_data_file)
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

def generate_sky_model_alms(fits_file,lmax=10):
    # http://healpy.readthedocs.org/en/latest/generated/healpy.sphtfunc.map2alm.html#healpy.sphtfunc.map2alm
    healmap = a.map.Map(fromfits=fits_file)
    as_pos = hp.sphtfunc.map2alm(healmap.map.map, lmax=lmax, pol=False)
    alms_pos = n.zeros([as_pos.shape[0],3])
    #print alms_pos.shape
    kk=0
    for ll in range(lmax+1):
        for mm in range(0,ll+1):
            alms_pos[kk] = n.array([ll,mm,as_pos[kk]])
            kk+=1
    #print alms_pos
    alms = n.zeros([(lmax+1)**2,3])
    kk=0
    for ll in range(lmax+1):
        for mm in range(-ll,ll+1):
            if mm<0:
                alm = alms_pos[n.where(n.logical_and(alms_pos[:,0]==ll, alms_pos[:,1]==-mm)),2]
                #print 'less',ll,mm,alm
                alms[kk] = n.array([ll,mm,n.conj(alm[0][0])])
            else:
                alm = alms_pos[n.where(n.logical_and(alms_pos[:,0]==ll, alms_pos[:,1]==mm)),2]
                #print 'greater ',ll,mm,alm 
                alms[kk] = n.array([ll,mm,alm[0][0]])
            kk+=1
    return alms 

def test_recover_alms(Q,N,a):
    # a is the alms from generate_sky_model_alms
    #print Q.shape, n.matrix(a).T.shape
    VV = n.array(n.matrix(Q)*n.matrix(a).T) + N 
    ahat,covar = uf.general_lstsq_fit_with_err(a,VV,Q,N,pseudo=5)
    err = n.sqrt(n.array(covar)*n.array(n.identity(covar.shape[0])))
    print "true      gs = {0}\nrecovered gs = {1}".format(a[0],ahat[0][0])
    #print "Error = ",err
    return a[0], ahat[0][0], err[0,0]

def compare_grids(lmax=2):
    beam_sigs = n.array([0.175,0.349,0.689,1.047]) #0.087,
    del_bls = n.array([4,6]) #,8,10,20
    param_grid = n.zeros([len(beam_sigs),len(del_bls)])#n.meshgrid(beam_sigs,del_bls)

    for ii,beam_sig in enumerate(beam_sigs):
        for jj,del_bl in enumerate(del_bls):
            Qstuff = n.load('./Q_matrices/grid_del_bl_{0:.2f}_num_bl_10_beam_sig_{1:.2f}_Q_max_l_10.npz'.format(del_bl,beam_sig))
            Q = Qstuff['Q']
            Q = Q[:,0:(lmax+1)**2]
            lms = Qstuff['lms']
            lms = lms[0:(lmax+1)**2,:]
            baselines = Qstuff['baselines']
            num_bl = len(baselines)
            print Q.shape, lms.shape, baselines.shape
            #N = total_noise_covar(0.1,num_bl,'./gsm_matrices/grid_del_bl_{0:.2f}_num_bl_10_beam_sig_{2:.2f}_gsm_max_l_10.npz'.format(del_bl,num_bl,beam_sig))
            N = (1.0**2)*n.identity(num_bl)
            alms = generate_sky_model_alms('/Users/mpresley/soft/gsm/haslam408_32.fits',lmax=lmax)
            gs_true, gs_recov, err = test_recover_alms(Q,N,alms[:,2])
            param_grid[ii,jj] = n.abs((gs_true-n.real(gs_recov))/gs_true)

    print param_grid.shape
    p.imshow(param_grid,interpolation='nearest',aspect='auto',extent=[0,len(del_bls),0,len(beam_sigs)]) #extent=[4,7,0.175,1.1]
    p.title('Log10 of Fractional difference btwn true and recovered global signal')
    p.yticks(range(len(beam_sigs)),beam_sigs)
    p.xticks(range(len(del_bls)),del_bls)
    p.ylabel('beam sigmas')
    p.xlabel('del baselines')
    p.colorbar()
    p.savefig('./figures/compare_grids.pdf')
    p.clf()

def monte_carlo(num):
    hasmap = a.map.Map(fromfits='/Users/mpresley/soft/gsm/haslam408_32.fits')
    gs_true = n.zeros(num)
    gs_recov = n.zeros(num)
    for ii in range(num):
        jiggled_map = haslam_extrap(hasdat=hasmap.map.map)
        jiggled_alms = hp.sphtfunc.map2alm(jiggled_map, lmax=10, pol=False)
        gs_true[ii], gs_recov[ii], _ = test_recover_alms(Q,N,jiggled_alms) #WRONG


if __name__=='__main__':
    Qstuff = n.load('./Q_matrices/grid_del_bl_4.00_num_bl_10_beam_sig_0.09_Q_max_l_10.npz')
    Q = Qstuff['Q']
    lms = Qstuff['lms']
    baselines = Qstuff['baselines']
    num_bl = len(baselines)#Q.shape[0]
    print 'bl ',num_bl
    N = total_noise_covar(0.1,num_bl,'./gsm_matrices/gsm_hybrid_del_bl_0.80_num_bl_7_lgbm_1.0_smbm_0.25.npz')
    # N = N*(4*n.pi/(3145728)**2)**2
    Q_N_plots(Q,N,lms,'grid')
    y = generate_sky_model_y(baselines,1.0,'/Users/mpresley/soft/gsm/haslam408_extrap_fq_0.1_32.fits')
    print y
    bl = n.sqrt(baselines[:,0]**2 + baselines[:,1]**2 + baselines[:,2]**2)
    p.scatter(bl,n.abs(y))
    p.savefig('./figures/visibilities_grid.pdf')
    # p.clf()
    # a = generate_sky_model_alms('./gsm_matrices/gsm_hybrid_del_bl_0.80_num_bl_7_lgbm_1.0_smbm_0.25.npz')
    # test_recover_alms(Q,N,a)

