import aipy as a, numpy as n, pylab as p
from scipy.stats import norm
import capo as C
import matplotlib as mpl
import healpy as hp
import useful_functions as uf
import global_sky_model as gsm
import sph_harm_coeffs as shc 
import Q_gsm_error_analysis as qgea 

"""
The monte carlo mpi code writes out a bunch of matrices. Each matrix has rows that 
are vectors of visibilities, so the number of columns is the number of baselines. 
Each row is a different set of visibilities, and each file has 10,000 rows.
"""

global data_loc; data_loc = '/Users/mpresley/Research/Research_Adrian_Aaron/nersc_data'#/matrices'
gsm_fits_file = '/Users/mpresley/soft/gsm/haslam408_extrap_fq_0.1_32.fits'#'/Users/mpresley/soft/gsm/data_50MHz_100MHz/gsm1001_32.fits'

def load_mc_data(mc_mat_loc):#del_bl=4.,num_bl=10,beam_sig=0.09,fq=0.05):
    #mc_mat_loc = '{0}/del_bl_{1:.2f}_num_bl_{2}_beam_sig_{3:.2f}'.format(mc_loc,del_bl,num_bl,beam_sig)
    for ii in range(1):
        stuff = n.load('{0}/mc_{1}.npz'.format(mc_mat_loc,ii))
        if ii==0:
            ys = stuff['matrix']
        else:
            ys = n.vstack(ys,stuff['matrix'])
    return ys 

def forground_hist(del_bl=8.,num_bl=10,beam_sig=0.09,fq=0.1):
    save_tag = 'grid_del_bl_{0:.2f}_num_bl_{1}_beam_sig_{2:.2f}_fq_{3:.2f}'.format(del_bl,num_bl,beam_sig,fq)
    print save_tag
    ys = load_mc_data('{0}/monte_carlo/{1}'.format(data_loc,save_tag))
    a00s = n.mean(n.abs(n.real(ys)),axis=1)

    _,bins,_ = p.hist(a00s,bins=36,normed=True)

    # plot best fit line
    mu,sigma = norm.fit(a00s)
    print "mu, sigma = ",mu,', ',sigma
    y_fit = mpl.mlab.normpdf(bins,mu,sigma)
    p.plot(bins, y_fit, 'r--', linewidth=2)

    p.xlabel('a_00')
    p.ylabel('Probability')
    p.title(save_tag)
    p.annotate('mu = {0:.2f}\nsigma = {1:.2f}'.format(mu,sigma), xy=(0.05, 0.95), xycoords='axes fraction')
    p.savefig('./figures/monte_carlo/true_sky_{0}.pdf'.format(save_tag))
    p.clf()



def construct_gs_hist(del_bl=8.,num_bl=10,beam_sig=0.09,fq=0.1):
    save_tag = 'grid_del_bl_{0:.2f}_num_bl_{1}_beam_sig_{2:.2f}_fq_{3:.2f}'.format(del_bl,num_bl,beam_sig,fq)
    ys = load_mc_data('{0}/monte_carlo/{1}'.format(data_loc,save_tag))
    
    alms_fg = qgea.generate_sky_model_alms(gsm_fits_file,lmax=3)
    alms_fg = alms_fg[:,2]

    baselines,Q,lms = qgea.load_Q_file(gh='grid',del_bl=del_bl,num_bl=num_bl,beam_sig=beam_sig,fq=fq,lmax=3)
    N = qgea.total_noise_covar(0.1,baselines.shape[0],'{0}/gsm_matrices/{1}.npz'.format(data_loc,save_tag))

    ahat00s = n.array([])
    for ii in xrange(ys.shape[0]):
        _,ahat,_ = qgea.test_recover_alms(ys[ii,:],Q,N,alms_fg,num_remov=None)
        ahat00s = n.append(n.real(ahat[0]),ahat00s)
    print ahat00s.shape
    _,bins,_ = p.hist(ahat00s,bins=36,normed=True)

    # plot best fit line
    mu,sigma = norm.fit(ahat00s)
    print "mu, sigma = ",mu,', ',sigma
    y_fit = mpl.mlab.normpdf(bins,mu,sigma)
    p.plot(bins, y_fit, 'r--', linewidth=2)

    p.xlabel('ahat_00')
    p.ylabel('Probability')
    p.title(save_tag)
    p.annotate('mu = {0:.2f}\nsigma = {1:.2f}'.format(mu,sigma), xy=(0.05, 0.95), xycoords='axes fraction')
    p.savefig('./figures/monte_carlo/{0}.pdf'.format(save_tag))
    p.clf()

if __name__=='__main__':
    for del_bl in (8.,10.,20.):
        for beam_sig in (0.09,0.17,0.35,0.69,1.05):
            forground_hist(del_bl=del_bl,beam_sig=beam_sig)
            construct_gs_hist(del_bl=del_bl,beam_sig=beam_sig)
