
#fqs = np.linspace(.14,.16, 50)
import matplotlib.pyplot as plt, numpy as np
import seaborn as sns
from w_opp import OppSolver
def freq_compare():
    sns.set_context(context='paper', font_scale=2)
    Afqs = np.linspace(.1, .2, 203)
    RES = []; RES2 = []; FQS = [np.array([Afqs[120]]), np.array([Afqs[140]])]
    RES3 = []; RES4 = [];
    C = ['m', 'c', 'b']
    bp = np.load('../64code/bandpass.npz')['bandpass']
    
    T1 = np.arange(2456681.35, 2456681.65, 0.002)
    DT = T1-np.mean(T1)
    for fqs in FQS:
    #fqs = Afqs[110:161]
    
        WS = OppSolver(fqs=fqs, bandpass=None, T1=T1)
        #res = WS.opp(bl1=(0,26), bl2=(0,26), return_series=True)
        #res2 = WS.opp(bl1=(0,26), bl2=(0,38), return_series=True)
        res = WS.opp(bl1=(0,103), bl2=(0,103), return_series=True, rephase=False)
        res2 = WS.opp(bl1=(0,103), bl2=(0,95), return_series=True, rephase=False)
        RES.append(res)
        RES2.append(res2)
        res3 = WS.opp(bl1=(0,103), bl2=(0,103), return_series=True, rephase=True)
        res4 = WS.opp(bl1=(0,103), bl2=(0,95), return_series=True, rephase=True)
        RES3.append(res3)
        RES4.append(res4)
        NORMS = [np.amax(np.abs(res)) for res in RES]
    #import IPython; IPython.embed()
    plt.figure()
    ax1 = plt.subplot(221)
    for i, fqs in enumerate(FQS):
        label = "%.2fGHz" % fqs[0]
        plt.plot(DT,RES[i].real/NORMS[i], c=C[i], label=label)
        plt.plot(DT,RES[i].imag/NORMS[i], '-.', c=C[i])
    plt.legend(loc=2)
    ax1.title.set_text('Without Rephasing')
    ax2 = plt.subplot(222)
    for i, fqs in enumerate(FQS):
        label = "%.2fGHz" % fqs[0]
        plt.plot(DT,RES3[i].real/NORMS[i], c=C[i], label=label)
        plt.plot(DT,RES3[i].imag/NORMS[i], '-.', c=C[i])
    plt.legend(loc=2)
    ax2.title.set_text('With Rephasing')
    ax3 = plt.subplot(223)
    for i, fqs in enumerate(FQS):
        label = "%.2fGHz" % fqs[0]
        plt.plot(DT,RES2[i].real/NORMS[i], c=C[i], label=label)
        plt.plot(DT,RES2[i].imag/NORMS[i], '-.', c=C[i]) 
        plt.plot(DT, np.abs(RES2[i])/NORMS[i], c=C[i], linewidth=5, alpha=0.5)
    plt.legend(loc=2) 
    plt.xlabel('Offset (Sidereal Days)')
    
    ax4 = plt.subplot(224)
    for i, fqs in enumerate(FQS):
        label = "%.2fGHz" % fqs[0]
        plt.plot(DT,RES4[i].real/NORMS[i], c=C[i], label=label)
        plt.plot(DT,RES4[i].imag/NORMS[i], '-.', c=C[i]) 
        plt.plot(DT, np.abs(RES4[i])/NORMS[i], c=C[i], linewidth=5, alpha=0.5)
    plt.legend(loc=2) 
    plt.xlabel('Offset (Sidereal Days)')
    plt.show()