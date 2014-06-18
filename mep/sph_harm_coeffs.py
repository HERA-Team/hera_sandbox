import aipy as a, numpy as n, pylab as p
import capo as C
import useful_functions as uf
from scipy import special
import matplotlib as mpl

def get_coeffs_lm(calfile,l,m,na = 32,freqs = n.array([.1,]),savefolderpath=None):
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

def get_Q(calfile,min_l,max_l,savefolderpath=None):
    for l in range(min_l,max_l+1):
        for m in range(-l,l+1):
            print l,m
            baselines,freqs,coeffs = get_coeffs_lm(calfile,l,m,na = 2,freqs = n.array([.1,]))
            if l==min_l and m==-min_l: Q = coeffs
            else: Q = n.hstack((Q,coeffs))
            if l==min_l and m==-min_l: lms = n.array([l,m])
            else: lms = n.vstack((lms,n.array([l,m])))
    print Q.shape
    if savefolderpath!=None: n.savez_compressed(savefolderpath+'{0}_Q_min_l_{1}_max_l_{2}'.format(calfile,min_l,max_l),Q=Q,baselines=baselines,lms=lms)
    return Q,baselines,lms

def plot_Q():#Q,lms):
    # n.savez_compressed('./coeff_data/basic_amp_aa_Q_min_l_0_max_l_8',

    Qstuff = n.load('./coeff_data/basic_amp_aa_Q_min_l_0_max_l_10.npz')
    Q = Qstuff['Q']
    lms = Qstuff['lms']
    p.scatter(lms[:,0],n.absolute(Q[0,:]),c=lms[:,1],cmap=mpl.cm.PiYG,s=50)
    p.yscale('log')
    p.ylim([10**-5,10**0.2])
    p.xlabel('l (color is m)')
    p.ylabel('Q for a baseline (1,1,0)')
    p.colorbar()
    p.show()

def combine_Q(Q1file,Q2file,newfile):
    Q1stuff = n.load(Q1file)
    Q2stuff = n.load(Q2file)
    Q = n.hstack((Q1stuff['Q'],Q2stuff['Q']))
    lms = n.vstack((Q1stuff['lms'],Q2stuff['lms']))
    baselines = n.vstack((Q1stuff['baselines'],Q2stuff['baselines']))
    n.savez_compressed(newfile,Q=Q,baselines=baselines,lms=lms)
    return Q, baselines, lms

if __name__=='__main__':
    calfile='basic_amp_aa'
    #na=8
    #baselines,freqs,coeffs = get_coeffs_lm(calfile,0,0,na=na,freqs=n.array([.1,]))
    #print coeffs
    #Q,baselines,lms = get_Q(calfile,9,10,savefolderpath='./coeff_data/')
    # Q, baselines, lms = combine_Q('./coeff_data/basic_amp_aa_Q_min_l_0_max_l_4.npz',
    #                             './coeff_data/basic_amp_aa_Q_min_l_5_max_l_8.npz',
    #                             './coeff_data/basic_amp_aa_Q_min_l_0_max_l_8')
    # Q, baselines, lms = combine_Q('./coeff_data/basic_amp_aa_Q_min_l_0_max_l_8.npz',
    #                             './coeff_data/basic_amp_aa_Q_min_l_9_max_l_10.npz',
    #                             './coeff_data/basic_amp_aa_Q_min_l_0_max_l_10')
    plot_Q()#Q,lms)
