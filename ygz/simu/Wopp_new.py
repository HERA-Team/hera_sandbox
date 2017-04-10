import aipy as a, numpy as np, capo as C, pylab as plt
from scipy import signal
import seaborn as sns
import w_opp
sns.set_context("paper", font_scale=2)


#bl1, bl2 = (0,103),(0,95)
#cuedict = {'26_23':0.125, '26_38': 0.042, '26_50': 0.083,'26_26':0., '50_57':0.122}
cuedict = {'26_26':0.,'26_38': 0.032557,'26_46': -0.034, '26_50':0.073557,'13_32':0.030557,'13_14':0.066557,'50_59':0.071557}
#MAXRES_EQUIV = 0.000150070063408  #peak of bl1=bl2
MAXRES_EQUIV = 1260.7792508                #peak of equivalent Opp
#MAXRES_EQUIV = 1.  #to compute MAXRES_EQUIV
#MAXRES_EQUIV = 2.51664842232e-05 #nside=128
#BLUENORM=0.18755
COMPARE = True
T0 = 2455700.5
T1 = np.arange(2455700.3,2455700.701,0.001)
WS = w_opp.OppSolver(fqs=np.array([.15]), T1=T1)
SERIES = []
for bls in (((0,103),(0,103)),((0,103),(0,95))):
	bl1, bl2 = bls
	if COMPARE:
	    fname = 'blout_'+str(bl1[0])+'_'+str(bl1[1])+'_'+str(bl2[0])+'_'+str(bl2[1])+'.npz'
	    print 'Reading', fname
	    file = np.load(fname)
	    TT, meancorr = file['arr_0'],file['arr_1']
	res = WS.opp(bl1=bl1,bl2=bl2, return_series=True)[:,0]
	SERIES.append(((T1,res),(TT,meancorr)))

norm=np.amax(np.abs(SERIES[0][0][1]))
sns.set(style="ticks", color_codes=True,font='DejaVu Serif', font_scale=2)
plt.rc('axes', linewidth=2)
f,axes = plt.subplots(2,1)
for i, obj in enumerate(SERIES):
	T1,res = SERIES[i][0]
	res /= norm
	axes[i].set_xlim([-0.2,0.2])
	axes[i].plot(T1-T0,res.real,'b',label='real')
	axes[i].plot(T1-T0,res.imag,'g',label='imag')
	axes[i].plot(T1-T0,np.abs(res),'r',alpha=0.5,linewidth=3,label='Theory')
	if COMPARE:
		TT, meancorr = SERIES[i][1]
		meancorr = meancorr/1.2875
		axes[i].plot(TT,meancorr.real,'b--')
		axes[i].plot(TT,meancorr.imag,'g--')
		axes[i].plot(TT,np.abs(meancorr),'r--',alpha=0.5,linewidth=3,label='Simulation')
	
	axes[i].get_yaxis().set_tick_params(direction='in')
	axes[i].get_xaxis().set_tick_params(direction='in')
	#plt.axvline(T1ac,color='k',alpha=0.5,linewidth=2)
axes[1].set_xlabel('dT (Sidereal Days)')
axes[0].legend(loc=2)
axes[0].set_title('Correlation Normalized to Equivalent Baseline Peak')
plt.setp(axes[0].get_xticklabels(), visible=False)
#plt.gcf().subplots_adjust(hspace=0.2)

plt.show()