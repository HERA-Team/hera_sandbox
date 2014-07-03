import aipy as a, numpy as n, pylab as p
import capo as C
import useful_functions as uf
import global_sky_model as gsm
from scipy import special
import matplotlib as mpl

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
            print ii,jj,[bx,by,bz,bb]
            baselines[ll] = [bx,by,bz,bb]
            ll+=1
    sorted_baselines = baselines[baselines[:,-1].argsort()] # sorts by last column
    sorted_baselines = sorted_baselines[:,0:3]
    if nb!=None: sorted_baselines = sorted_baselines[0:nb,:]
    return sorted_baselines

def get_coeffs_lm_fewer_baselines(aa,baselines,l,m,freqs = n.array([.1,]),savefolderpath=None,beamsig=None):
    """
    This function calculates the coefficients in front of the lm spherical
    harmonic for the antenna array described in the calfile. The coefficients 
    are determined by the integral of A(l,m)Y(l,m)exp[-i2pi(ul+vm)]dldm.

    However, this function only loops over the baselines of the antennas with
    respect to the origin. 
    """
    im = a.img.Img(size=200, res=.5) #make an image of the sky to get sky coords
    tx,ty,tz = im.get_top(center=(200,200)) #get coords of the zenith?
    valid = n.logical_not(tx.mask)
    tx,ty,tz = tx.flatten(),ty.flatten(),tz.flatten()
    theta = n.arctan(ty/tx) # using math convention of theta=[0,2pi], phi=[0,pi]
    phi = n.arccos(n.sqrt(1-tx*tx-ty*ty))
    #n.set_printoptions(threshold='nan')
    #print theta,phi
    #quit()

    #beam response for an antenna pointing at (tx,ty,tz) with a polarization in x direction
    #amp = A(theta) in notes
    if beamsig==None:
        amp = aa[0].bm_response((tx,ty,tz),pol='x')**2
    else:
        amp = uf.gaussian(beamsig,n.zeros_like(theta),phi) 
    na = len(aa.ants) # number of antennas
    #coefficient array: rows baselines; cols frequencies
    coeffs = n.zeros([(na*na-na)/2,len(freqs)],dtype=n.complex)
    # compute normalization for spherical harmonics
    Ynorm = special.sph_harm(0,0,0,0)
    # loop through all baselines
    for jj in range(baselines.shape[0]):
        bx,by,bz = baselines[jj]
        kk=0
        for fq in freqs: #loop over frequencies
            phs = n.exp(-2j*n.pi*fq * (bx*tx+by*ty+bz*tz)) #fringe pattern
            Y = n.array(special.sph_harm(m,l,theta,phi))/Ynorm #using math convention of theta=[0,2pi], phi=[0,pi]
            Y.shape = phs.shape = amp.shape = im.uv.shape
            amp = n.where(valid, amp, 0)
            phs = n.where(valid, phs, 0)
            Y = n.where(valid, Y, 0) 
            # n.set_printoptions(threshold='nan')
            # print Y
            # quit() 
            dc_response = n.sum(amp*Y*phs)/n.sum(amp) #this integrates the amplitude * fringe pattern; units mK?
            print 'dc = ',dc_response
            jy_response = n.real(dc_response * 100 / C.pspec.jy2T(fq)) # jy2T converts flux density in jansky to temp in mK
            print '\t',fq, dc_response, jy_response
            coeffs[jj,kk] = dc_response
            kk+=1
    if savefolderpath!=None: n.savez_compressed('{0}{1}_data_l_{1}_m_{2}'.format(savefolderpath,calfile,l,m),baselines=baselines,frequencies=freqs,coeffs=coeffs)
    return baselines,freqs,coeffs


def get_coeffs_lm(calfile,l,m,freqs = n.array([.1,]),savefolderpath=None):
    """
    This function calculates the coefficients in front of the lm spherical
    harmonic for the antenna array described in the calfile. The coefficients 
    are determined by the integral of A(l,m)Y(l,m)exp[-i2pi(ul+vm)]dldm.
    """
    aa = a.cal.get_aa(calfile, n.array([.150])) #get antenna array
    im = a.img.Img(size=200, res=.5) #make an image of the sky to get sky coords
    tx,ty,tz = im.get_top(center=(200,200)) #get coords of the zenith?
    valid = n.logical_not(tx.mask)
    tx,ty,tz = tx.flatten(),ty.flatten(),tz.flatten()
    theta = n.arctan(ty/tx) # using math convention of theta=[0,2pi], phi=[0,pi]
    phi = n.arccos(n.sqrt(1-tx*tx-ty*ty))
    #n.set_printoptions(threshold='nan')
    #print theta,phi
    #quit()

    #beam response for an antenna pointing at (tx,ty,tz) with a polarization in x direction
    #amp = A(theta) in notes
    amp = aa[0].bm_response((tx,ty,tz),pol='x')**2 
    na = len(aa.ants) # number of antennas
    #coefficient array: rows baselines; cols frequencies
    coeffs = n.zeros([(na*na-na)/2,len(freqs)],dtype=n.complex)
    # compute normalization for spherical harmonics
    Ynorm = special.sph_harm(0,0,0,0)
    # loop through all baselines
    baselines = n.zeros([(na*na-na)/2,3])
    ll=0
    for ii in n.arange(na):
        for jj in n.arange(ii):
            bx,by,bz = aa.get_baseline(ii,jj,'z') 
            #the baseline for antennas ii and jj 
            print ii,jj,[bx,by,bz]
            baselines[ll] = [bx,by,bz]
            kk=0
            for fq in freqs: #loop over frequencies
                phs = n.exp(-2j*n.pi*fq * (bx*tx+by*ty+bz*tz)) #fringe pattern
                Y = n.array(special.sph_harm(m,l,theta,phi))/Ynorm #using math convention of theta=[0,2pi], phi=[0,pi]
                Y.shape = phs.shape = amp.shape = im.uv.shape
                amp = n.where(valid, amp, 0)
                phs = n.where(valid, phs, 0)                
                Y = n.where(valid, Y, 0) 
                # n.set_printoptions(threshold='nan')
                # print Y
                # quit() 
                dc_response = n.sum(amp*Y*phs)/n.sum(amp) #this integrates the amplitude * fringe pattern; units mK?
                print 'dc = ',dc_response
                jy_response = n.real(dc_response * 100 / C.pspec.jy2T(fq)) # jy2T converts flux density in jansky to temp in mK
                print '\t',fq, dc_response, jy_response
                coeffs[ll,kk] = dc_response
                kk+=1
            ll+=1
    if savefolderpath!=None: n.savez_compressed('{0}{1}_data_l_{1}_m_{2}'.format(savefolderpath,calfile,l,m),baselines=baselines,frequencies=freqs,coeffs=coeffs)
    return baselines,freqs,coeffs

def get_coeffs_lm_baselines_from_origin(calfile,l,m,freqs = n.array([.1,]),savefolderpath=None):
    """
    This function calculates the coefficients in front of the lm spherical
    harmonic for the antenna array described in the calfile. The coefficients 
    are determined by the integral of A(l,m)Y(l,m)exp[-i2pi(ul+vm)]dldm.

    However, this function only loops over the baselines of the antennas with
    respect to the origin. 
    """
    aa = a.cal.get_aa(calfile, n.array([.150])) #get antenna array
    im = a.img.Img(size=200, res=.5) #make an image of the sky to get sky coords
    tx,ty,tz = im.get_top(center=(200,200)) #get coords of the zenith?
    valid = n.logical_not(tx.mask)
    tx,ty,tz = tx.flatten(),ty.flatten(),tz.flatten()
    theta = n.arctan(ty/tx) # using math convention of theta=[0,2pi], phi=[0,pi]
    phi = n.arccos(n.sqrt(1-tx*tx-ty*ty))
    #n.set_printoptions(threshold='nan')
    #print theta,phi
    #quit()

    #beam response for an antenna pointing at (tx,ty,tz) with a polarization in x direction
    #amp = A(theta) in notes
    amp = aa[0].bm_response((tx,ty,tz),pol='x')**2 
    na = len(aa.ants) # number of antennas
    #coefficient array: rows baselines; cols frequencies
    coeffs = n.zeros([(na*na-na)/2,len(freqs)],dtype=n.complex)
    # compute normalization for spherical harmonics
    Ynorm = special.sph_harm(0,0,0,0)
    # loop through all baselines
    baselines = n.zeros([(na*na-na)/2,3])
    ll=0
    for jj in range(1,na):
        bx,by,bz = aa.get_baseline(0,jj,'z') 
        #the baseline for antennas 0 and jj 
        print 0,jj,[bx,by,bz]
        baselines[ll] = [bx,by,bz]
        kk=0
        for fq in freqs: #loop over frequencies
            phs = n.exp(-2j*n.pi*fq * (bx*tx+by*ty+bz*tz)) #fringe pattern
            Y = n.array(special.sph_harm(m,l,theta,phi))/Ynorm #using math convention of theta=[0,2pi], phi=[0,pi]
            Y.shape = phs.shape = amp.shape = im.uv.shape
            amp = n.where(valid, amp, 0)
            phs = n.where(valid, phs, 0)                
            Y = n.where(valid, Y, 0) 
            # n.set_printoptions(threshold='nan')
            # print Y
            # quit() 
            dc_response = n.sum(amp*Y*phs)/n.sum(amp) #this integrates the amplitude * fringe pattern; units mK?
            print 'dc = ',dc_response
            jy_response = n.real(dc_response * 100 / C.pspec.jy2T(fq)) # jy2T converts flux density in jansky to temp in mK
            print '\t',fq, dc_response, jy_response
            coeffs[ll,kk] = dc_response
            kk+=1
        ll+=1
    if savefolderpath!=None: n.savez_compressed('{0}{1}_data_l_{1}_m_{2}'.format(savefolderpath,calfile,l,m),baselines=baselines,frequencies=freqs,coeffs=coeffs)
    return baselines,freqs,coeffs

def get_Q(calfile,min_l,max_l,mvals=None,nb=6,savefolderpath=None,beamsig=None):
    """
    This function creates a Q matrix for the antenna array in the calfile.
    Note that although a 'full' Q matrix would vary l from 0 to max_l, this
    function gives the option of starting from a higher l. This enables you 
    to calculate the Q matrix in chunks since the computation takes a while
    (about 33 sec per element in the Q matrix).
    """
    aa = a.cal.get_aa(calfile, n.array([.150])) #get antenna array
    baselines = get_baselines(aa,nb=nb)
    if mvals==None: fullm=True
    else: fullm=False
    for l in range(min_l,max_l+1):
        if fullm: mvals = range(-l,l+1)
        for m in mvals:
            print l,m
            baselines,freqs,coeffs = get_coeffs_lm_fewer_baselines(aa,baselines,l,m,freqs=n.array([.1,]),beamsig=beamsig)
            print 'got coeffs'
            if l==min_l and m==mvals[0]: Q = coeffs
            else: Q = n.hstack((Q,coeffs))
            if l==min_l and m==mvals[0]: lms = n.array([l,m])
            else: lms = n.vstack((lms,n.array([l,m])))
    print Q.shape
    if savefolderpath!=None: n.savez_compressed(savefolderpath+'{0}_Q_min_l_{1}_max_l_{2}'.format(calfile,min_l,max_l),Q=Q,baselines=baselines,lms=lms)
    return Q,baselines,lms

def plot_Q(Q,lms,save_tag=None):
    # n.savez_compressed('./coeff_data/basic_amp_aa_Q_min_l_0_max_l_8',

    p.scatter(lms[:,0],n.absolute(Q[0,:]),c=lms[:,1],cmap=mpl.cm.PiYG,s=50)
    p.yscale('log')
    p.ylim([10**-5,10**0.2])
    p.xlabel('l (color is m)')
    p.ylabel('First row of Q matrix')
    p.colorbar()
    p.savefig('./figures/{0}_10_Q_matrix_elements.pdf'.format(save_tag))
    p.clf()
    #p.show()

def combine_Q(Q1file,Q2file,newfile):
    Q1stuff = n.load(Q1file)
    Q2stuff = n.load(Q2file)
    Q = n.hstack((Q1stuff['Q'],Q2stuff['Q']))
    lms = n.vstack((Q1stuff['lms'],Q2stuff['lms']))
    #baselines = n.vstack((Q1stuff['baselines'],Q2stuff['baselines']))
    baselines = Q1stuff['baselines']
    n.savez_compressed(newfile,Q=Q,baselines=baselines,lms=lms)
    return Q, baselines, lms

def combine_Q_baselines(Q1file,Q2file,newfile):
    Q1stuff = n.load(Q1file)
    Q2stuff = n.load(Q2file)
    Q = n.vstack((Q1stuff['Q'],Q2stuff['Q']))
    baselines = n.vstack((Q1stuff['baselines'],Q2stuff['baselines']))
    n.savez_compressed(newfile,Q=Q,baselines=baselines,lms=lms)
    return Q, baselines, lms

def test_recover_gs(Q,baselines,lms,n_sig=5):
    print lms[:,0].shape
    a = get_a_from_gsm(lms[:,0]) # this is the sky, i.e. the x in y=Qx+n 
    a.shape = (a.shape[0],1)
    print Q.shape
    print a.shape
    VV = n.array(n.matrix(Q)*n.matrix(a))+n.random.normal(loc=0.0,scale=n_sig,size=[Q.shape[0],1])*n.exp(2*n.pi*1j*n.random.rand())
    print a.shape
    print VV.shape
    ahat,covar = uf.general_lstsq_fit_with_err(a,VV,Q,(n_sig**2)*n.identity(Q.shape[0]),pseudo=True)
    print covar
    err = n.sqrt(n.array(covar)*n.array(n.identity(covar.shape[0])))
    print err.shape
    print "true      gs = {0}\nrecovered gs = {1}".format(a[0],ahat[0])
    print "Error = ",err
    return a[0], ahat[0], err[0,0]

def test_recover_gs_vary_n(Q,baselines,lms):
    """
    This function runs many tests of the linear regression, varying the 
    amount of noise introduced into the data. I hard-coded the global 
    signal strength to be 1 so that it is easy to compare the magnitude
    of the recovered and true global signals and the amplitude of the noise.
    """
    for jj in n.arange(100):
        gs_diff = n.zeros(20,dtype=complex)
        errors = n.zeros(20)
        n_sigs = n.logspace(-3,1,num=20)
        print n_sigs
        for ii,n_sig in enumerate(n_sigs):
            print ii
            gs_true,gs_recov,err = test_recover_gs(Q,baselines,lms,n_sig=n_sig)
            print gs_true
            print gs_recov
            print gs_true.shape
            gs_diff[ii] = gs_recov[0] - gs_true[0] 
            errors[ii] = err
        p.scatter(n_sigs,n.absolute(gs_diff))
        p.scatter(n_sigs,errors,color="red")
        p.xscale('log')
        p.yscale('log')
        p.xlim(1e-4,1e2)
    p.xlabel('Amplitude of noise relative to global signal\n(I.e. true global signal amplitude is 1)')
    #p.ylabel('Recovered global signal (true gs = 1)')
    p.ylabel('Difference between true and recovered global signal')
    #p.show()
    p.savefig('./figures/circle10_Q_pinv_gs_diff_vs_n.pdf')
    p.clf()

def get_a_from_gsm(l):
    print l.shape
    C = n.where(l<=8, n.exp(-1.450*l+0.1003*l*l), 0.7666*(l**(-2.365))) # gsm power_spectrum from http://arxiv.org/pdf/1404.2596v2.pdf pg. 13
    print C.shape
    a = n.zeros_like(C,dtype='complex')
    for ii in range(len(l)):
        a[ii] = n.random.normal(loc=0.0,scale=C[ii]/2) + 1j*n.random.normal(loc=0.0,scale=C[ii]/2)
    return a

def info_matrix(Q,N,lms,save_tag=None):
    Q = n.matrix(Q); N = n.matrix(N)
    info = Q.H*N.I*Q
    #print info
    
    p.scatter(lms[:,0],n.diag(info),c=lms[:,1],cmap=mpl.cm.PiYG,s=50)
    p.xlabel('l (color is m)')
    p.ylabel('diagonal elements of info matrix')
    # p.yscale('log')
    # p.ylim([10**-5,10])
    p.colorbar()
    p.savefig('./figures/{0}_info_matrix_diagonal.pdf'.format(save_tag))
    p.clf()
    
    first_row = n.array(n.absolute(info[0,:]))
    p.scatter(lms[:,0],first_row,c=lms[:,1],cmap=mpl.cm.PiYG,s=50)
    p.xlabel('l (color is m)')
    p.ylabel('Abs val of first row of info matrix')
    # p.yscale('log')
    # p.ylim([10**-5,10])
    p.colorbar()
    p.savefig('./figures/{0}_info_matrix_first_row.pdf'.format(save_tag))
    p.clf()



    eig_vals,eig_vecs = n.linalg.eig(info)
    p.scatter(range(len(eig_vals)),eig_vals,s=50)
    p.xlabel('index of eigenvalue array (I dont think they are in any particular order)')
    p.ylabel('eigenvalue')
    # p.yscale('log')
    # p.ylim([10**-10,10**-2])
    p.savefig('./figures/{0}_info_matrix_eig_vals.pdf'.format(save_tag))
    p.clf()
    
    foo = n.array(n.absolute(eig_vecs[:,0]))
    p.scatter(lms[:,0],foo,c=lms[:,1],cmap=mpl.cm.PiYG,s=50)
    p.xlabel('l (color is m)')
    p.ylabel('coefficient of a_l,m for the largest eigenvector')
    # p.yscale('log')
    # p.ylim([10**-5,10])
    p.savefig('./figures/{0}_info_matrix_eig_vec.pdf'.format(save_tag))
    print eig_vals.shape
    print 'eig_vals = ',eig_vals[0]
    print 'eig_vecs = ',eig_vecs[:,0]

def window_function_matrix(Q,N,lms,save_tag=None):
    ls = n.arange(max(lms[:,0])+1)
    l_locs = ls*ls

    M = n.matrix(n.zeros_like(Q))
    Q = n.matrix(Q); N = n.matrix(N)
    Ninv = N.I
    info = Q.H*Ninv*Q
    
    p.imshow(n.log(n.absolute(info)))
    p.title('Info Matrix')
    p.xticks(l_locs,ls)
    p.yticks(l_locs,ls)
    p.xlabel('l')
    p.ylabel('l')
    p.savefig('./figures/{0}_info_matrix_ufpseudo.pdf'.format(save_tag))
    #p.show()
    p.clf()

    # for ii in range(M.shape[0]):
    #     M[ii,ii] = 1/info[ii,ii]
    #M = n.linalg.pinv(info)
    M = uf.pseudo_inverse(info,num_remov=1)
    #n.set_printoptions(threshold='nan')
    #print M*info
    #BB = n.absolute(n.absolute(M*info)-n.identity(M.shape[0]))
    #print n.amax(BB)
    W = M*info
    #print M[0:4,0:4]
    #print info[0:4,0:4]
    #print W[0:4,0:4]
    p.imshow(n.log(n.absolute(W)))
    p.title('Window Function Matrix')
    p.xticks(l_locs,ls)
    p.yticks(l_locs,ls)
    p.xlabel('l')
    p.ylabel('l')
    p.savefig('./figures/{0}_W_matrix_ufpseudo.pdf'.format(save_tag))
    #p.show()
    p.clf()

    foo = n.array(n.absolute(W[0,:]))
    p.scatter(lms[:,0],foo,c=lms[:,1],cmap=mpl.cm.PiYG,s=50)
    #p.yscale('log')
    #p.ylim([10**-3,10**0.2])
    p.xlabel('l (color is m)')
    p.ylabel('first row of Window Function Matrix')
    p.colorbar()
    p.savefig('./figures/{0}_W_pinv_matrix_elements.pdf'.format(save_tag))
    #p.show()
    p.clf()


def fringe_pattern_plots(baselines,lms):
    im = a.img.Img(size=200, res=.5) #make an image of the sky to get sky coords
    tx,ty,tz = im.get_top(center=(200,200)) #get coords of the zenith?
    valid = n.logical_not(tx.mask)
    tx,ty,tz = tx.flatten(),ty.flatten(),tz.flatten()
    theta = n.arctan(ty/tx) # using math convention of theta=[0,2pi], phi=[0,pi]
    phi = n.arccos(n.sqrt(1-tx*tx-ty*ty))
    Ynorm = special.sph_harm(0,0,0,0)
    fq = 0.1

    # for jj in range(lms.shape[0]):
    #     l,m = lms[jj]
    #     Y = n.array(special.sph_harm(m,l,theta,phi))/Ynorm #using math convention of theta=[0,2pi], phi=[0,pi]
    #     Y.shape = im.uv.shape
    #     p.imshow(n.real(Y))
    #     p.title('Spherical Harmonic l,m = {0}, {1}'.format(l,m))
    #     p.xlabel('tx')
    #     p.ylabel('ty')
    #     p.savefig('./figures/fringe_patterns/Y_l_{0}_m_{1}.pdf'.format(l,m))
    #     #p.show()
    #     p.clf()

    for jj in range(baselines.shape[0]):
        bx,by,bz = baselines[jj,:]
        phs = n.exp(-2j*n.pi*fq * (bx*tx+by*ty+bz*tz)) #fringe pattern
        phs.shape = im.uv.shape
        p.imshow(n.real(phs))
        p.title('Fringe Pattern bx,by,bz = {0:.2f}, {1:.2f}, {2:.2f}'.format(bx,by,bz))
        p.xlabel('tx')
        p.ylabel('ty')
        p.savefig('./figures/fringe_patterns/phs_bx_{0:.2f}_by_{1:.2f}_bz_{2:.2f}.pdf'.format(bx,by,bz))
        #p.show()
        p.clf()

    # for beamsig in (5,10,15,20,25,30):
    #     beamsig_rad = beamsig*n.pi/180.
    #     amp = uf.gaussian(beamsig_rad,theta,phi) 
    #     amp.shape = im.uv.shape
    #     p.imshow(n.absolute(amp))
    #     p.title('Gaussian Beam for sigma = {0}'.format(beamsig))
    #     p.xlabel('tx')
    #     p.ylabel('ty')
    #     p.savefig('./figures/fringe_patterns/bm_sig_{0}.pdf'.format(beamsig))
    #     #p.show()
    #     p.clf()


if __name__=='__main__':
    #baselines,freqs,coeffs = get_coeffs_lm(calfile,0,0,freqs=n.array([.1,]))
    #print coeffs
    #Q,baselines,lms = shc.get_Q('basic_amp_aa_long',3000,3001,mvals=(0,1500,3000),savefolderpath='./coeff_data/long/')
    #Q,baselines,lms = shc.get_Q('basic_amp_aa_long',3500,3501,mvals=(0,1750,3500),savefolderpath='./coeff_data/long/')
    
    #Q,baselines,lms = shc.get_Q('basic_amp_aa_circle_gauss',0,3,savefolderpath='./coeff_data/circle_13_gauss_30deg/',beamsig=0.524)
    #Q,baselines,lms = shc.get_Q('basic_amp_aa_circle_gauss',4,5,savefolderpath='./coeff_data/circle_13_gauss_30deg/',beamsig=0.524)
    #Q,baselines,lms = shc.get_Q('basic_amp_aa_circle_gauss',6,7,savefolderpath='./coeff_data/circle_13_gauss_30deg/',beamsig=0.524)

    #Q,baselines,lms = shc.get_Q('basic_amp_aa_circle_gauss',8,9,savefolderpath='./coeff_data/circle_13_gauss_15deg/',beamsig=0.262)
    #Q,baselines,lms = shc.get_Q('basic_amp_aa_circle_gauss',10,11,savefolderpath='./coeff_data/circle_13_gauss_15deg/',beamsig=0.262)
    #Q,baselines,lms = shc.get_Q('basic_amp_aa_circle_gauss',12,13,savefolderpath='./coeff_data/circle_13_gauss_15deg/',beamsig=0.262)



    keyword = 'hybrid_grid_Q_max'
    Q, baselines, lms = combine_Q('./Q_matrices/hybrid_grid_1_Q_max_l_15.npz',
                                './Q_matrices/hybrid_grid_2_Q_max_l_15.npz',
                                './Q_matrices/hybrid_grid_12_Q_max_l_15')
    Q, baselines, lms = combine_Q('./Q_matrices/hybrid_grid_12_Q_max_l_15.npz',
                                './Q_matrices/hybrid_grid_3_Q_max_l_15.npz',
                                './Q_matrices/hybrid_grid_123_Q_max_l_15')
    Q, baselines, lms = combine_Q('./Q_matrices/hybrid_grid_123_Q_max_l_15.npz',
                                './Q_matrices/hybrid_grid_4_Q_max_l_15.npz',
                                './Q_matrices/hybrid_grid_Q_max_l_15')

    #calfile='basic_amp_aa_circle'
    # Q, baselines, lms = combine_Q('./coeff_data/{0}/basic_amp_aa_{1}_Q_min_l_0_max_l_3.npz'.format(keyword,keyword),
    #                             './coeff_data/{0}/basic_amp_aa_{1}_Q_min_l_4_max_l_5.npz'.format(keyword,keyword),
    #                             './coeff_data/{0}/basic_amp_aa_{1}_Q_min_l_0_max_l_5'.format(keyword,keyword))
    # Q, baselines, lms = combine_Q('./coeff_data/{0}/basic_amp_aa_{1}_Q_min_l_0_max_l_5.npz'.format(keyword,keyword),
    #                             './coeff_data/{0}/basic_amp_aa_{1}_Q_min_l_6_max_l_7.npz'.format(keyword,keyword),
    #                             './coeff_data/{0}/basic_amp_aa_{1}_Q_min_l_0_max_l_7'.format(keyword,keyword))
    # Q, baselines, lms = combine_Q('./coeff_data/{0}/basic_amp_aa_{1}_Q_min_l_0_max_l_7.npz'.format(keyword,keyword),
    #                             './coeff_data/{0}/basic_amp_aa_{1}_Q_min_l_8_max_l_9.npz'.format(keyword,keyword),
    #                             './coeff_data/{0}/basic_amp_aa_{1}_Q_min_l_0_max_l_9'.format(keyword,keyword))
    # Q, baselines, lms = combine_Q('./coeff_data/{0}/basic_amp_aa_{1}_Q_min_l_0_max_l_9.npz'.format(keyword,keyword),
    #                             './coeff_data/{0}/basic_amp_aa_{1}_Q_min_l_10_max_l_11.npz'.format(keyword,keyword),
    #                             './coeff_data/{0}/basic_amp_aa_{1}_Q_min_l_0_max_l_11'.format(keyword,keyword))
    # Q, baselines, lms = combine_Q('./coeff_data/{0}/basic_amp_aa_{1}_Q_min_l_0_max_l_11.npz'.format(keyword,keyword),
    #                             './coeff_data/{0}/basic_amp_aa_{1}_Q_min_l_12_max_l_13.npz'.format(keyword,keyword),
    #                             './coeff_data/{0}/basic_amp_aa_{1}_Q_min_l_0_max_l_13'.format(keyword,keyword))
    #Qstuff = n.load('./coeff_data/{0}/basic_amp_aa_{1}_Q_min_l_0_max_l_7.npz'.format(keyword,keyword))
    Qstuff = n.load('./Q_matrices/hybrid_grid_Q_max_l_15.npz')
    Q = Qstuff['Q']
    lms = Qstuff['lms']
    baselines = Qstuff['baselines']

    # baselines = baselines[0,:]
    # Q = Q[0,:]
    # Q.shape = n.array([1,Q.shape[0]])
    # baselines.shape = n.array([1,baselines.shape[0]])
    # print Q.shape
    # print baselines.shape

    keyword = keyword

    #print baselines
    #fringe_pattern_plots(baselines,lms)
    plot_Q(Q,lms,save_tag=keyword)
    #aa = a.cal.get_aa(calfile, n.array([.10]))
    #amp = aa[0].bm_response((500,100,1000),pol='x')**2 
    #print amp
    # Nfg = gsm.gsm_noise_covar(baselines,aa,savepath='./coeff_data/{0}/gsm_noise_covar'.format(keyword))
    # p.imshow(n.log(n.absolute(Nfg)))
    # p.show()

    N = (1.0**2)*n.identity(Q.shape[0])
    window_function_matrix(Q,N,lms,save_tag=keyword)
    info_matrix(Q,N,lms,save_tag=keyword)
    #test_recover_gs_vary_n(Q,baselines,lms)
    #test_recover_gs(Q,baselines,lms,n_sig=.1)
