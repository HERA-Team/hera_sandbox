
#fqs = np.linspace(.14,.16, 50)
import pylab as plt
import seaborn as sns
def freq_compare():
    sns.set_context(context='paper', font_scale=2)
    Afqs = np.linspace(.1, .2, 203)
    RES = []; RES2 = []; FQS = [np.array([Afqs[120]]), np.array([Afqs[140]])]
    C = ['m', 'c', 'b']
    bp = np.load('../64code/bandpass.npz')['bandpass']
    
    T1 = np.arange(2456681.35, 2456681.65, 0.002)
    DT = T1-np.mean(T1)
    for fqs in FQS:
    #fqs = Afqs[110:161]
    
        WS = OppSolver(fqs=fqs, bandpass=None, T1=T1)
        res = WS.opp(bl1=(0,103), bl2=(0,103), return_series=True)
        res2 = WS.opp(bl1=(0,103), bl2=(0,95), return_series=True)
        RES.append(res)
        RES2.append(res2)
        NORMS = [np.amax(np.abs(res)) for res in RES]
    #import IPython; IPython.embed()
    plt.figure()
    plt.subplot(211)
    for i, fqs in enumerate(FQS):
        label = "%.2fGHz" % fqs[0]
        plt.plot(DT,RES[i].real/NORMS[i], c=C[i], label=label)
        plt.plot(DT,RES[i].imag/NORMS[i], '-.', c=C[i])
    plt.legend()
    plt.subplot(212)
    for i, fqs in enumerate(FQS):
        label = "%.2fGHz" % fqs[0]
        plt.plot(DT,RES2[i].real/NORMS[i], c=C[i], label=label)
        plt.plot(DT,RES2[i].imag/NORMS[i], '-.', c=C[i]) 
        plt.plot(DT, np.abs(RES2[i])/NORMS[i], c=C[i], linewidth=5, alpha=0.5)
    plt.legend() 
    plt.xlabel('Offset (Sidereal Days)')
    plt.show()