#! /usr/bin/env python

import aipy as a, numpy as n, pylab as p
from scipy.stats import norm
import matplotlib as mpl
import useful_functions as uf

global data_loc; data_loc = "/global/homes/a/acliu/globalSig/fq_120_150_testCase"

def load_Q_file(savekey,fq,lmax=-1):
    Q_file = '{0}/Q_matrices/Q_{1}_fq_{2:.3f}.npz'.format(data_loc,savekey,fq)
    
    Qstuff = n.load(Q_file)
    Q = Qstuff['Q']
    lms = Qstuff['lms']
    baselines = Qstuff['baselines']
    if lmax != -1:
        Q = Q[:,0:(lmax+1)**2]
        lms = lms[0:(lmax+1)**2,:]
    Qstuff.close()
    return baselines,Q,lms

def total_noise_covar(savekey,eps,fq,Nfg_type):
    if Nfg_type == "gsm":
        fg_file = '{0}/gsm_matrices/gsm_{1}_fq_{2:.3f}.npz'.format(data_loc,savekey,fq)
    if Nfg_type == "improved":
        fg_file = '{0}/gsm_matrices/improvedNfg_{1}_fq_{2:.3f}.npz'.format(data_loc,savekey2,fq)

    Nfg_file = n.load(fg_file)
    Nfg = Nfg_file['matrix']
    Nfg_file.close()
    Ninst = (eps**2)*n.identity(Nfg.shape[0])
    Ntot = Nfg + Ninst 
    return Ntot

def return_MQdagNinv(Q,N,mode,savekey,fq,Nfg_type,eps=10**-4,num_remov=None):
    #Q = n.matrix(Q); N = n.matrix(N)
    #Ninv = uf.pseudo_inverse(N,num_remov=None) # XXX want to remove dynamically
    Ninv = n.linalg.inv(N)
    info = n.dot(n.conjugate(Q.T),n.dot(Ninv,Q))
    #info = n.dot(Q.H,n.dot(Ninv,Q))
    eigens = n.linalg.eigvalsh(info).real
    diags = n.diagonal(info).real

    plotTitleBase = '{0}_fq_{1:.3f}'.format(savekey,fq)
    f, ax = p.subplots(figsize=(8,6), dpi=400)
    ax.plot(eigens)
    ax.set_ylabel("Information matrix eigenvalues", fontsize=20, fontname="Times New Roman")
    ax.set_xlabel("Eigenvalue number", fontsize=20, fontname="Times New Roman")
    ax.set_title('Info_scree_\nN_fg_{1}_\n{0}'.format(plotTitleBase,Nfg_type))
    f.savefig('{0}/plots/{1}/Info_scree_Nfg_{3}_{2}.pdf'.format(data_loc,savekey,plotTitleBase,Nfg_type))
                                                                                                                                                                                                                    
    f, ax = p.subplots(figsize=(8,6), dpi=400)
    ax.plot(diags)
    ax.set_ylabel("Information matrix diagonals", fontsize=20, fontname="Times New Roman")
    ax.set_xlabel("Element number", fontsize=20, fontname="Times New Roman")
    ax.set_title('Info_diagonal_\nN_fg_{1}_\n{0}'.format(plotTitleBase,Nfg_type))
    f.savefig('{0}/plots/{1}/Info_diagonal_Nfg_{3}_{2}.pdf'.format(data_loc,savekey,plotTitleBase,Nfg_type))

    if mode == "pseudo":
        M = uf.pseudo_inverse(info,num_remov=num_remov)
        MQN = n.dot(M,n.dot(n.conjugate(Q.T),Ninv))
        MQN = n.array(MQN) 
    elif mode == "regularization":
        info += eps**2 * n.identity(info.shape[0])
        MQN = n.dot(n.linalg.inv(info),n.dot(n.conjugate(Q.T),Ninv))
    elif mode == "diagonal":
        diag_info = n.diag(n.diag(info))
        MQN = n.dot(n.linalg.inv(diag_info),n.dot(n.conjugate(Q.T),Ninv))

    WindMatrix = n.dot(MQN,Q)
    WindNorm = WindMatrix[0,0]
    WindMatrixRow = WindMatrix[0] / WindNorm
    MQN_firstRowNormed = MQN[0] / WindNorm

    return MQN_firstRowNormed, WindMatrixRow


def analyze(savekey,fqs,mode,Nfg_type,Ntot_eps=10**-4,info_eps=10**-4):
    T0s_master = []
    for fq in fqs:
        currentFolder = '{0}/MCs/{1}_fq_{2:.3f}'.format(data_loc,savekey,fq)
        print "Looking in...",currentFolder
        ys = n.load('{0}/mc_{1}_fq_{2:.3f}_allMCs.npz'.format(currentFolder,savekey,fq))['matrix']
        print "....loaded",ys.shape[0],"MCs, each of which has info about",ys.shape[1],"baselines."
        print "Now loading Q matrices..."
        baselines, Q, lms = load_Q_file(savekey,fq,-1)
        print "...loaded.  Onto loading the N_fg matrices..."
        Ntot = total_noise_covar(savekey,Ntot_eps,fq,Nfg_type)
        MQN_firstRow, WindMatrixRow = return_MQdagNinv(Q,Ntot,mode,savekey,fq,Nfg_type,info_eps,num_remov=None)
        ahat00s = n.einsum('j,ij',MQN_firstRow,ys)
        T0s_master.append(ahat00s.real/n.sqrt(4.*n.pi))
    T0s_master = n.array(T0s_master) # fqs x MCs
    numMCs = T0s_master.shape[1]
    T0s_bias = n.mean(T0s_master,axis=1)
    T0s_covar = n.zeros((len(fqs),len(fqs)))
    for i,spec in enumerate(T0s_master.T):
        T0s_covar += n.outer(spec.real,spec.real)
    T0s_covar /= numMCs
    T0s_covar -= n.outer(T0s_bias.real,T0s_bias.real)

    # Save the data
    n.savez_compressed('{0}/MCs/T0_MCs_Nfg_type_{3}_{2}_{1}'.format(data_loc,savekey,mode,Nfg_type),T0s=T0s_master)
    n.save('{0}/MCs/T0_bias_Nfg_type_{3}_{2}_{1}.npy'.format(data_loc,savekey,mode,Nfg_type),T0s_bias)
    n.save('{0}/MCs/T0_covar_Nfg_type_{3}_{2}_{1}.npy'.format(data_loc,savekey,mode,Nfg_type),T0s_covar)

    # Plot window functions
    WindMatrixRow = n.real(n.array(WindMatrixRow))
    f, ax = p.subplots(figsize=(8,6), dpi=400)
    im = ax.scatter(lms[:,0],WindMatrixRow,c=lms[:,1],cmap=mpl.cm.PiYG,s=50)
    ax.set_xlabel('l (color is m)')
    ax.set_ylabel('Monopole window function')
    ax.set_title(savekey)
    f.colorbar(im)
    f.savefig('{0}/plots/{1}/Window_Nfg_type_{3}_{2}_{1}.pdf'.format(data_loc,savekey,mode,Nfg_type))


    # Plot histograms
    for T0s,fq in zip(T0s_master,fqs):
        plotTitleBase = '{0}_fq_{1:.3f}'.format(savekey,fq)
        f, ax = p.subplots(figsize=(8,6), dpi=400)
        _,bins,_ = p.hist(T0s,bins=36,normed=True)
        mu,sigma = norm.fit(T0s)
        y_fit = mpl.mlab.normpdf(bins,mu,sigma)
        ax.plot(bins, y_fit, 'r--', linewidth=2)
        ax.set_xlabel(r'$\widehat{T}_0$')
        ax.set_ylabel(r'$P(\widehat{T}_0)$')
        ax.set_title(plotTitleBase)
        p.annotate('mu = {0:.2f}\nsigma = {1:.2f}'.format(mu,sigma), xy=(0.05, 0.95), xycoords='axes fraction')
        f.savefig('{0}/plots/{1}/T0_histogram_Nfg_type_{4}_{3}_{2}.pdf'.format(data_loc,savekey,plotTitleBase,mode,Nfg_type))

    # Plot spectra
    errorBars = n.sqrt(n.diag(T0s_covar))
    f, ax = p.subplots(figsize=(8,6), dpi=400)
    ax.errorbar(fqs, T0s_bias, yerr=errorBars, fmt='o',ecolor='Black')
    ax.set_title(savekey)
    ax.set_xlim([fqs[0]-0.001,fqs[-1]+0.001])
    ax.set_xlabel('Frequency [GHz]')
    ax.set_ylabel(r'Recovered $T_0$')
    f.savefig('{0}/plots/{1}/T0_spectrum_Nfg_type_{3}_{2}_{1}.pdf'.format(data_loc,savekey,mode,Nfg_type))

    return None

if __name__=='__main__':
    beam_sigs = [1.57] #(np.pi/18,np.pi/6,5*np.pi/18,7*np.pi/18)
    sqGridSideLens = [12] #(4,8,12,16)
    variableBeams = [0] #(0,1)
    lowerFreq = 120.
    upperFreq = 150. #150.
    freqSpace = 1.
    fqs = n.arange(lowerFreq,upperFreq+freqSpace,freqSpace)
    fqs /= 1000. # Convert from MHz to GHz
    lowerFreq /= 1000.
    mode = "diagonal"
    Ntot_eps = 10**-4
    info_eps = 10**-4

    for beam_sig in beam_sigs:
        del_bl = 1/(2*n.pi*beam_sig*lowerFreq)
        for sqGridSideLen in sqGridSideLens:
            for variableBeam in variableBeams:
                if variableBeam == 0:
                    savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_fixedWidth_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)
                elif variableBeam == 1:
                    savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_lambdaBeam_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)
                #savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)
                #savekey2 = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_fixedWidth_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)

                analyze(savekey,fqs,mode,"gsm",Ntot_eps,info_eps)
                #analyze(savekey,fqs,mode,"improved",Ntot_eps,info_eps)





