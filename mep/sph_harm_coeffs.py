import aipy as a, numpy as n, pylab as p
import capo as C
import useful_functions as uf
from scipy import special

def get_coeffs_lm(calfile,l,m,na = 32,freqs = n.arange(.1,)):
    """
    This function calculates the coefficients in front of the lm spherical
    harmonic for the antenna array described in the calfile. The coefficients 
    are determined by the integral of A(l,m)Y(l,m)exp[-i2pi(ul+vm)]dldm.

    Not working yet. Doesn't give same results for l,m=0,0 as previous get_coeff code.
    """
    aa = a.cal.get_aa(calfile, n.array([.150])) #get antenna array
    im = a.img.Img(size=200, res=.5) #make an image of the sky to get sky coords
    tx,ty,tz = im.get_top(center=(200,200)) #get coords of the zenith?
    valid = n.logical_not(tx.mask)
    tx,ty,tz = tx.flatten(),ty.flatten(),tz.flatten()
    print tx,ty
    theta = n.arcsin(ty/tx) # using math convention of theta=[0,2pi], phi=[0,pi]
    phi = n.arccos(n.sqrt(1-tx*tx-ty*ty))

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
                Y = special.sph_harm(m,l,theta,phi)/Ynorm #using math convention of theta=[0,2pi], phi=[0,pi]
                print Y
                Y.shape = phs.shape = amp.shape = im.uv.shape 
                amp = n.where(valid, amp, 0)
                phs = n.where(valid, phs, 0)
                Y = n.where(valid, Y, 0)

                dc_response = n.sum(amp*Y*phs)/n.sum(amp) #this integrates the amplitude * fringe pattern; units mK?
                jy_response = n.real(dc_response * 100 / C.pspec.jy2T(fq)) # jy2T converts flux density in jansky to temp in mK
                print '\t',fq, dc_response, jy_response
                coeffs[ll,kk] = dc_response
                kk+=1
            ll+=1
    n.savez_compressed('./coeff_data/{0}_data_l_{1}_m_{2}'.format(calfile,l,m),baselines=baselines,frequencies=freqs,coeffs=coeffs)
    return baselines,freqs,coeffs

def get_Q(calfile,max_l,savefolderpath=None):
    for l in range(max_l+1):
        for m in range(-max_l-1,max_l+1):
            print l,m
            baselines,freqs,coeffs = get_coeffs_lm(calfile,l,m,na = 32,freqs = n.arange(.1,))
            if l==m==0: Q = coeffs
            else: Q = n.hstack((Q,coeffs))
    print Q.shape
    if savefolderpath!=None: n.savez_compressed(savefolderpath+'{0}_Q_max_l_{1}'.format(calfile,max_l),Q=Q,baselines=baselines)
    return Q

if __name__=='__main__':
    calfile='basic_amp_aa'
    na=8
    baselines,freqs,coeffs = get_coeffs_lm(calfile,0,0,na=na,freqs=n.array([.1,]))
    print coeffs


