import aipy as a
import sys
import capo
from glob import glob
import pylab as plt
import numpy as np

def flag_avg(data,flags):
    avg = (data*flags).sum(axis=0)/flags.sum(axis=0)
    return avg

def tuplist2str(tuplist):
    string = ''
    for tup in tuplist:
        string += str(tup[0])+'_'+str(tup[1])+','
    string = string[0:-1]
    return string

dir = '/Users/jaguirre/Documents/PAPER/GlobalSignal/'
files = sorted(glob(dir+'zen.2456680.[2,3,4,5]*.uvcRREO'))
#[dir+'zen.2456680.17382.xx.uvcRREO',dir+'zen.2456680.18078.xx.uvcRREO']
#sys.argv[1:]

anttuple = [(64,72),(3,61),(49,73)]
antstrings = tuplist2str(anttuple)

#antstrings = '3_41,18_76,7_94,3_2,41_90,55_7' #,1_2,2_3'
#anttuple = [(3,41),(18,76),(7,94),(2,3),(41,90),(7,55)]

polstr = 'xx'

print "getting data"
tinfo,wfall,flags = capo.arp.get_dict_of_uv_data(files,antstr=antstrings,polstr=polstr)
print "done getting data"

# The waterfall is wfall[(1,2)]['xx'] 
#%%
w = {}
iflags = {}
for anttup in anttuple:
    w[anttup] = wfall[anttup]['xx']
    ifl = np.ones(flags[anttup]['xx'].shape)
    ifl[flags[anttup]['xx']] = 0
    iflags[anttup] = ifl

#%%

spec = {}
for anttup in anttuple:
    plt.figure()
    plt.subplot(121)
    plt.imshow(np.log10(np.abs(w[anttup])),aspect='auto')
    plt.colorbar()
    plt.subplot(122)
    plt.imshow(np.angle(w[anttup]),aspect='auto')
    spec[anttup] = flag_avg(w[anttup],iflags[anttup])

#%%
MHz = np.linspace(100,200,num=len(spec[spec.keys()[0]]))
plt.figure()
for anttup in anttuple:
    plt.subplot(211)
    plt.plot(MHz,spec[anttup].imag,'--',label=str(anttup))
    plt.subplot(212)
    plt.plot(MHz,spec[anttup].real,label=str(anttup))
plt.subplot(211)
plt.legend()
plt.subplot(212)
plt.legend()