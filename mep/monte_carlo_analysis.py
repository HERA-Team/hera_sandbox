import aipy as a, numpy as n, pylab as p
import Q_gsm_error_analysis as qgea
import useful_functions as uf
from scipy.stats import norm
import matplotlib as mpl
import os
import healpy as hp 
from load_data import *

"""
The monte carlo mpi code writes out a bunch of matrices. Each matrix has rows that 
are vectors of visibilities, so the number of columns is the number of baselines. 
Each row is a different set of visibilities, and each file has 10,000 rows.
"""

def return_MQdagNinv(Q,N,num_remov=None):
    Q = n.matrix(Q); N = n.matrix(N)
    Ninv = uf.pseudo_inverse(N,num_remov=None) # XXX want to remove dynamically
    info = n.dot(Q.H,n.dot(Ninv,Q))
    M = uf.pseudo_inverse(info,num_remov=num_remov)
    MQN = n.dot(M,n.dot(Q.H,Ninv))
    #print n.diag(N) 
    return MQN 

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
    save_tag = 'grid_del_bl_{0:.2f}_num_bl_{1}_beam_sig_{2:.2f}_fq_{3:.3f}'.format(del_bl,num_bl,beam_sig,fq)
    save_tag_mc = 'grid_del_bl_{0:.2f}_num_bl_{1}_beam_sig_{2:.2f}_fq_{3}'.format(del_bl,num_bl,beam_sig,fq)
    ys = load_mc_data('{0}/monte_carlo/{1}'.format(data_loc,save_tag_mc))
    print 'ys ',ys.shape
    
    alms_fg = qgea.generate_sky_model_alms(gsm_fits_file,lmax=3)
    alms_fg = alms_fg[:,2]

    baselines,Q,lms = load_Q_file(gh='grid',del_bl=del_bl,num_bl=num_bl,beam_sig=beam_sig,fq=fq,lmax=3)
    N = total_noise_covar(0.1,baselines.shape[0],'{0}/gsm_matrices/gsm_{1}.npz'.format(data_loc,save_tag))
    MQN = return_MQdagNinv(Q,N,num_remov=None)
    print MQN
    ahat00s = n.array([])
    for ii in xrange(ys.shape[1]):
        #_,ahat,_ = qgea.test_recover_alms(ys[:,ii],Q,N,alms_fg,num_remov=None)
        ahat = uf.vdot(MQN,ys[:,ii])
        ahat00s = n.append(n.real(ahat[0]),ahat00s)
    #print ahat00s
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
    p.annotate('mu = {0:.2f}\nsigma = {1:.2f}'.format(mu,sigma), xy=(0.05, 0.5), xycoords='axes fraction')
    p.savefig('./figures/monte_carlo/{0}.pdf'.format(save_tag))
    p.clf()

def construct_covar(del_bl=8.,num_bl=10,beam_sig=0.09,savea00=True):
    save_tag_base = 'grid_del_bl_{0:.2f}_num_bl_{1}_beam_sig_{2:.2f}'.format(del_bl,num_bl,beam_sig)
    print save_tag_base
    ahat_mat = n.array([]); fqs = n.array([])
    for mc_fold in os.listdir('{0}/monte_carlo'.format(data_loc)):
        if save_tag_base in mc_fold:
            ahat00s = n.load('{0}/monte_carlo/{1}/ahat00s.npz'.format(data_loc,mc_fold))
            ahat00s = ahat00s['ahat00s']
            fq = float(mc_fold.split('_')[-1])
            ahat00s.shape = (1,ahat00s.shape[0])
            if ahat_mat.shape==(0,): ahat_mat = ahat00s
            else: ahat_mat = n.vstack((ahat_mat,ahat00s))
            fqs = n.append(fqs,fq)
    print ahat_mat.shape # num_fq x num_mc
    if savea00: n.savez_compressed('{0}/monte_carlo/spectrum_data/{1}_ahat00s'.format(data_loc,save_tag_base),ahat00s=ahat_mat,fqs=fqs)
    ahat00_avgs = n.mean(ahat_mat,axis=1) 
    covar = n.zeros((ahat_mat.shape[0],ahat_mat.shape[0]))
    for ii in range(ahat_mat.shape[0]):
        for jj in range(ii+1):
            covar[ii,jj] = covar[jj,ii] = n.mean(ahat_mat[ii,:]*ahat_mat[jj,:]) - ahat00_avgs[ii]*ahat00_avgs[jj]
    return covar, fqs 

def plot_covar(del_bl=8.,beam_sig=0.09,save_covar=True):
    covar,fqs = construct_covar(del_bl=del_bl,beam_sig=beam_sig)
    #print covar 
    p.pcolor(fqs,fqs,covar)
    p.colorbar()
    p.savefig('{0}/mc_spec_figs/{1}_covar.pdf'.format(fig_loc,save_tag_base))
    if save_covar: n.savez_compressed('{0}/monte_carlo/spectrum_data/{1}_covar'.format(data_loc,save_tag_base),covar=covar,fqs=fqs)

def plot_spectrum(del_bl=8.,num_bl=10,beam_sig=0.09,savea00=True,lmax=3):
    save_tag_base = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_beam_sig_{2:.2f}'.format(del_bl,num_bl,beam_sig)
    print save_tag_base

    save_tag_base_old = 'grid_del_bl_{0:.2f}_num_bl_{1}_beam_sig_{2:.2f}'.format(del_bl,num_bl,beam_sig)

    mu = n.array([]); sigma = n.array([]) 
    fqs = n.array([]); nums = n.array([])

    for mc_fold in os.listdir('{0}/monte_carlo'.format(data_loc)):
        if save_tag_base_old in mc_fold:
            fq = float(mc_fold.split('_')[-1])
    #for fq in (0.050,0.064,0.076,0.088,0.090):
            mc_fold = '{0}_fq_{1}'.format(save_tag_base,fq)
            ys = load_mc_data('{0}/monte_carlo'.format(data_loc,del_bl=del_bl,num_bl=num_bl,beam_sig=beam_sig,fq=fq))
            baselines,Q,lms = load_Q_file(gh='grid',del_bl=del_bl,num_bl=num_bl,beam_sig=beam_sig,fq=fq,lmax=lmax)
            #N = total_noise_covar(0.0,del_bl=del_bl,num_bl=num_bl,beam_sig=beam_sig,fq=fq)
            N = n.eye(baselines.shape[0])
            MQN = return_MQdagNinv(Q,N,num_remov=None)
            #print 'MQN \n',MQN 
            ahat00s = n.array([])
            for ii in xrange(ys.shape[1]):
                #print MQN.shape
                #print 'ys shape ',ys.shape
                ahat = uf.vdot(MQN,ys[:,ii])
                #print ahat[0]
                ahat00s = n.append(n.real(ahat[0]),ahat00s)
            mu0,sigma0 = norm.fit(ahat00s)
            mu = n.append(mu,mu0); sigma = n.append(sigma,sigma0)
            fqs = n.append(fqs,fq); nums = n.append(nums,ys.shape[1])
            print fq
            #if fq==0.06: print ahat00s 
            if savea00: n.savez_compressed('{0}/monte_carlo/{1}/ahat00s'.format(data_loc,mc_fold),ahat00s=ahat00s)
    print fqs
    print mu
    print sigma
    p.errorbar(fqs, mu, yerr=sigma, fmt='o',ecolor='Black')#,c=nums,cmap=mpl.cm.copper_r)
    #p.scatter(fqs,mu,color='g')
    p.title('Mean and sigma of recovered global signal for\n{0}'.format(save_tag_base))
    p.xlim([0.045,0.100])
    #p.ylim([-1500.,0.0])
    p.xlabel('Freq (GHz)')
    p.ylabel('Recovered a00')
    p.savefig('{0}/mc_spec_figs/{1}_lmax_{2}.pdf'.format(fig_loc,save_tag_base,lmax))
    p.clf()

def plot_haslam_spectrum(del_bl=8.,num_bl=10,beam_sig=0.09,savea00=True,lmax=3):
    save_tag_base = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_beam_sig_{2:.2f}'.format(del_bl,num_bl,beam_sig)
    print save_tag_base

    save_tag_base_old = 'grid_del_bl_{0:.2f}_num_bl_{1}_beam_sig_{2:.2f}'.format(del_bl,num_bl,beam_sig)

    mu = n.array([]); sigma = n.array([]) 
    fqs = n.array([]); nums = n.array([])

    for mc_fold in os.listdir('{0}/monte_carlo'.format(data_loc)):
        if save_tag_base_old in mc_fold:
            fq = float(mc_fold.split('_')[-1])
    #for fq in (0.050,0.064,0.076,0.088,0.090):
            mc_fold = '{0}_fq_{1}'.format(save_tag_base,fq)
            ys = load_mc_data('{0}/monte_carlo_haslam'.format(data_loc,del_bl=del_bl,num_bl=num_bl,beam_sig=beam_sig,fq=fq))
            baselines,Q,lms = load_Q_file(gh='grid',del_bl=del_bl,num_bl=num_bl,beam_sig=beam_sig,fq=fq,lmax=lmax)
            #N = total_noise_covar(0.0,del_bl=del_bl,num_bl=num_bl,beam_sig=beam_sig,fq=fq)
            N = n.eye(baselines.shape[0])
            MQN = return_MQdagNinv(Q,N,num_remov=None)
            #print 'MQN \n',MQN 
            ahat00s = n.array([])
            for ii in xrange(ys.shape[1]):
                #print MQN.shape
                #print 'ys shape ',ys.shape
                ahat = uf.vdot(MQN,ys[:,ii])
                #print ahat[0]
                ahat00s = n.append(n.real(ahat[0]),ahat00s)
            mu0,sigma0 = norm.fit(ahat00s)
            mu = n.append(mu,mu0); sigma = n.append(sigma,sigma0)
            fqs = n.append(fqs,fq); nums = n.append(nums,ys.shape[1])
            print fq
            print mu0 
            #if fq==0.06: print ahat00s 
            if savea00: n.savez_compressed('{0}/monte_carlo_haslam/{1}/ahat00s'.format(data_loc,mc_fold),ahat00s=ahat00s)
    print fqs
    print mu
    print sigma
    p.errorbar(fqs, mu, yerr=sigma, fmt='o',ecolor='Black')#,c=nums,cmap=mpl.cm.copper_r)
    #p.scatter(fqs,mu,color='g')
    p.title('Mean and sigma of recovered global signal for\n{0}'.format(save_tag_base))
    p.xlim([0.045,0.100])
    #p.ylim([-1500.,0.0])
    p.xlabel('Freq (GHz)')
    p.ylabel('Recovered a00')
    p.savefig('{0}/mc_spec_figs/{1}_haslam.pdf'.format(fig_loc,save_tag_base))
    p.clf()

if __name__=='__main__':
    for del_bl in (4.,):#6.,8.):
        for beam_sig in (0.09,):#0.17,0.35,0.69,1.05):
            num_bl = 10
            # save_tag_base = 'grid_del_bl_{0:.2f}_num_bl_{1}_beam_sig_{2:.2f}'.format(del_bl,num_bl,beam_sig)
            plot_spectrum(del_bl=del_bl,num_bl=num_bl,beam_sig=beam_sig,savea00=False,lmax=3)
            # for fq in (0.05,0.06,0.09,0.074):
            #     save_tag_base = 'W_grid_del_bl_{0:.2f}_sqGridSideLen_{1}_beam_sig_{2:.2f}_fq_{3}'.format(del_bl,num_bl,beam_sig,fq)
            #     baselines,Q,lms = load_Q_file(gh='grid',del_bl=del_bl,num_bl=num_bl,beam_sig=beam_sig,fq=fq,lmax=3)
            #     p.imshow(n.real(Q))
            #     p.colorbar()
            #     p.title('Q')
            #     p.savefig('{0}/Q_{1}.pdf'.format(fig_loc,save_tag_base))
            #     p.clf()
            #     #N = total_noise_covar(0.0,del_bl=del_bl,num_bl=num_bl,beam_sig=beam_sig,fq=fq)
            #     N = n.eye(baselines.shape[0])
            #     qgea.window_fn_matrix(Q,N,num_remov=None,save_tag=save_tag_base,lms=lms)
            #     construct_gs_hist(del_bl=del_bl,num_bl=num_bl,beam_sig=beam_sig,fq=fq)
            