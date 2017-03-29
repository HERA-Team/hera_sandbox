
#fqs = np.linspace(.14,.16, 50)
import matplotlib.pyplot as plt, numpy as np
import seaborn as sns
from w_opp import OppSolver
from matplotlib.font_manager import FontProperties
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
        res = WS.opp(bl1=(0,103), bl2=(0,103), return_series=True, rephase=0)
        res2 = WS.opp(bl1=(0,103), bl2=(0,95), return_series=True, rephase=0)
        RES.append(res)
        RES2.append(res2)
        res3 = WS.opp(bl1=(0,103), bl2=(0,103), return_series=True, rephase=104.936)
        res4 = WS.opp(bl1=(0,103), bl2=(0,95), return_series=True, rephase=104.936)
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

def freq_compare_v2():
    Afqs = np.linspace(.1, .2, 203)
    RES = []; RES2 = []; FQS = [np.array([Afqs[120]]), np.array([Afqs[140]])]
    RES3 = []; RES4 = [];
    
    bp = np.load('../64code/bandpass.npz')['bandpass']
    
    T1 = np.arange(2456681.35, 2456681.65, 0.002)
    DT = T1-np.mean(T1)
    for fqs in FQS:
    #fqs = Afqs[110:161]
    
        WS = OppSolver(fqs=fqs, bandpass=None, T1=T1)
        #res = WS.opp(bl1=(0,26), bl2=(0,26), return_series=True)
        #res2 = WS.opp(bl1=(0,26), bl2=(0,38), return_series=True)
        res = WS.opp(bl1=(0,103), bl2=(0,103), return_series=True, rephase=0)
        res2 = WS.opp(bl1=(0,103), bl2=(0,95), return_series=True, rephase=0)
        RES.append(res)
        RES2.append(res2)
        res3 = WS.opp(bl1=(0,103), bl2=(0,103), return_series=True, rephase=104.936)
        res4 = WS.opp(bl1=(0,103), bl2=(0,95), return_series=True, rephase=104.936)
        RES3.append(res3)
        RES4.append(res4)
        
    #import IPython; IPython.embed()
    return RES,RES2,RES3,RES4,FQS,DT

def plot_freq_compare(RES,RES2,RES3,RES4,FQS,DT):
    NORMS = [np.amax(np.abs(res)) for res in RES]
    #font = FontProperties()
    #font.set_weight('bold')
    #sns.set_context(context='paper')
    sns.set(style="ticks", color_codes=False,font='DejaVu Serif', font_scale=1.5)
    plt.rc('axes', linewidth=2)
    C = ['c', 'm']

    fig = plt.figure()
    ax1 = plt.subplot(411)
    for i, fqs in enumerate(FQS):
        label = "%.2fGHz" % fqs[0]
        plt.plot(DT,RES[i].real/NORMS[i], c=C[i], label=label, linewidth=2)
    ax1.set_ylim([-1.1,1.1])
        #plt.plot(DT,RES[i].imag/NORMS[i], '-.', c=C[i])
        #plt.plot(DT, np.abs(RES[i])/NORMS[i], c=C[i], linewidth=8, alpha=0.5)
    plt.legend(loc=2)
    ax1.set_ylabel('Without Rephasing')
    ax2 = plt.subplot(413)
    for i, fqs in enumerate(FQS):
        label = "%.2fGHz" % fqs[0]
        plt.plot(DT,RES3[i].real/NORMS[i], c=C[i], label=label, linewidth=2)
        #plt.plot(DT,RES3[i].imag/NORMS[i], '-.', c=C[i])
        #plt.plot(DT, np.abs(RES3[i])/NORMS[i], c=C[i], linewidth=8, alpha=0.5)
    ax2.set_ylim([-1.1,1.1])
    plt.legend(loc=2)
    ax2.set_ylabel('With Rephasing')
    ax3 = plt.subplot(412)
    for i, fqs in enumerate(FQS):
        label = "%.2fGHz" % fqs[0]
        plt.plot(DT,RES2[i].real/NORMS[i], c=C[i], label=label, linewidth=2)
        #plt.plot(DT,RES2[i].imag/NORMS[i], '-.', c=C[i]) 
        #plt.plot(DT, np.abs(RES2[i])/NORMS[i], c=C[i], linewidth=8, alpha=0.5)
    plt.legend(loc=2) 
    #plt.xlabel('Offset (Sidereal Days)')
    
    ax4 = plt.subplot(414)
    for i, fqs in enumerate(FQS):
        label = "%.2fGHz" % fqs[0]
        plt.plot(DT,RES4[i].real/NORMS[i], c=C[i], label=label, linewidth=2)
        #plt.plot(DT,RES4[i].imag/NORMS[i], '-.', c=C[i]) 
        #plt.plot(DT, np.abs(RES4[i])/NORMS[i], c=C[i], linewidth=8, alpha=0.5)

    for ax in [ax1,ax2,ax3,ax4]:
        if ax is not ax4:
            plt.setp(ax.get_xticklabels(), visible=False)
        start, end = ax.get_ylim()
        start+=0.1; end-=0.1
        print start, end
        plt.grid()
        ax.yaxis.set_ticks(np.linspace(start, end, 3))
    fig.subplots_adjust(hspace=-0.05)
    plt.legend(loc=2) 
    plt.xlabel('Offset (Sidereal Days)')
    plt.show()
