import pyuvdata as uv
import hera_pspec as ps
import numpy as np
import pylab as plt
from glob import glob
import capo as C
from capo import oqe
from itertools import cycle
from random import shuffle

# List of filenames of the data to load
dfiles_odd = glob('/Users/josh/Desktop/PSA64/odd/sep0,1/lst.*3.*.uvGLI')
dfiles_even = glob('/Users/josh/Desktop/PSA64/even/sep0,1/lst.*2.*.uvGLI')

class DataSet(oqe.DataSet):
    """Extention of oqe.DataSet to include some covariance regularization."""

    def iC(self, k, t=None, rcond=1e-12):
        """Regularize covariance before inverting."""
        assert(t is None)
        if not self._iC.has_key(k):
            C = self.C(k)
            try:
                C = C + np.identity(len(C))*np.trace(C)#*float(opts.mode_num)
            except(ValueError):
                print "Tried to use this mode_num: ", opts.mode_num
                print "using instead 0 regularization"
                C = C
            U, S, V = np.linalg.svd(C.conj())
            iS = 1./S
            self.set_iC({k: np.einsum('ij,j,jk', V.T, iS, U.T)})
        return self._iC[k]

    def set_data(self, dsets, wgts=None, conj=None):
        """Set data inside of object, also computes average over bls."""
        if type(dsets.values()[0]) == dict:
            dsets, wgts = self.flatten_data(dsets), self.flatten_data(wgts)
            self.x, self.w = {}, {}
        avgx = []
        avgC = []
        for k in dsets:
            self.x[k] = dsets[k].T
            try: self.w[k] = wgts[k].T
            except(TypeError): self.w[k] = np.ones_like(self.x[k])
            print k[1]
            try:
                if conj[k[1]]: self.x[k] = np.conj(self.x[k])
            except(TypeError, KeyError): pass
            try:
                avgx.append(self.x[k])
                avgC.append(oqe.cov(self.x[k], self.w[k]))
            except(TypeError, KeyError): pass
        self.avgx = np.average(avgx, axis=0)
        self.avgC = np.average(avgC, axis=0)


def make_PS(keys, dsv, dsn, grouping=False):
    """Use OQE formalism to generate power spectrum.                                                                                       
    Output weighted and identity weightings.                                                                                               
    """
    newkeys = keys
    # identity and covariance case dataset is the same                                                                                 
    dsIv, dsCv = dsv, dsv
    dsIn, dsCn = dsn, dsn
    pCvs, pIvs = [], []
    pCns, pIns = [], []
    for k, key1 in enumerate(newkeys):
        if k == 1 and len(newkeys) == 2:
            # NGPS = 1 (skip 'odd' with 'even' if we already did 'even' with 'odd')
            continue
        for key2 in newkeys[k:]:
            # Check for even/odd multiply and that it's the same bsl
            if key1[2] == key2[2] or key1[:2] != key2[:2]:
                print 'Skipped.'
                continue
            else:
                scalar = 1.
                FCv = dsCv.get_F(key1, key2, cov_flagging=False)
                FIv = dsIv.get_F(key1, key2, use_cov=False, cov_flagging=False)
                qCv = dsCv.q_hat(key1, key2, cov_flagging=False)
                qIv = dsIv.q_hat(key1, key2, use_cov=False, cov_flagging=False)
                MCv, WCv = dsCv.get_MW(FCv, mode='I')
                MIv, WIv = dsIv.get_MW(FIv, mode='I')
                pCv = dsCv.p_hat(MCv, qCv, scalar=scalar)
                pIv = dsIv.p_hat(MIv, qIv, scalar=scalar)
                pCvs.append(pCv)
                pIvs.append(pIv)

                FCn = dsCn.get_F(key1, key2, cov_flagging=False)
                FIn = dsIn.get_F(key1, key2, use_cov=False, cov_flagging=False)
                qCn = dsCn.q_hat(key1, key2, cov_flagging=False)
                qIn = dsIn.q_hat(key1, key2, use_cov=False, cov_flagging=False)
                MCn, WCn = dsCn.get_MW(FCn, mode='I')
                MIn, WIn = dsIn.get_MW(FIn, mode='I')
                pCn = dsCv.p_hat(MCn, qCn, scalar=scalar)
                pIn = dsIv.p_hat(MIn, qIn, scalar=scalar)
                pCns.append(pCn)
                pIns.append(pIn)
    return np.array(pIvs), np.array(pIns)

new_files_odd = []
new_files_even = []
lsts_odd = []
lsts_ODD = []
#### PAPER LST binned files dont have equal time samples, find files with 21 LST samples
for i in dfiles_odd:
    _d = uv.UVData()
    _d.read_miriad(i)
    if len(np.unique(_d.time_array)) == 21:
        new_files_odd.append(i)
        lsts_ODD.append(np.unique(_d.lst_array))
        lsts_odd.append(np.mean(np.unique(_d.lst_array)))

lsts_ODD = np.sort(np.array(lsts_ODD).reshape(-1))

lsts_even = []
lsts_EVEN = []
for i in dfiles_even:
    _d = uv.UVData()
    _d.read_miriad(i)
    if len(np.unique(_d.time_array)) == 21:
        new_files_even.append(i)
        lsts_even.append(np.mean(np.unique(_d.lst_array)))
        lsts_EVEN.append(np.unique(_d.lst_array))
lsts_EVEN = np.sort(np.array(lsts_EVEN).reshape(-1))


new_files_even = np.array(new_files_even)
new_files_odd = np.array(new_files_odd)

# Not really 'aligning', just getting close enough
align_even = list(new_files_even[np.argsort(lsts_even)])
align_odd = list(new_files_odd[np.argsort(lsts_odd)])
len_even = len(align_even)
align_odd = align_odd[:len_even]
print 'Odd/Even: ',len(dfiles_odd),len(dfiles_even)
print 'Pruned Odd/Even: ',len(align_odd),len(align_even)

# Load into UVData objects
d = []
et_ct = 0 # equal time count
eo_grp = 0
interlaced_files = []
bad_files = []

### Find the difference between valid even/odds for terminating
eo_diff = np.abs(len(dfiles_odd)-len(dfiles_even))
print eo_diff

_d = uv.UVData()
_d.read_miriad(align_odd)
_d.select(freq_chans=(range(95,116)),polarizations=1)
d.append(_d)
_d = uv.UVData()
_d.read_miriad(align_even)
_d.select(freq_chans=(range(95,116)),polarizations=1)
d.append(_d)

print 'Number of antennae',np.shape(d[0].get_antpairs())
# Specify which baselines to include
bls = d[0].get_antpairs()
permute_bls = [(a,b) for a in bls for b in bls]
cross_bls = []
for a,b in permute_bls:
    if a != b or a[::-1] == b:
        continue
    else:
        cross_bls.append((a,b))
print 'Shape of cross_bls: ',np.shape(cross_bls)

###### NEW HERA PSPEC OQE ######
w = [None]*2
ds = ps.PSpecData(dsets=d, wgts=w)
newHERApspec, pairs = ds.pspec(cross_bls, input_data_weight='identity', norm='I', verbose=True)
print 'Shape of New HERA OQE pspec: ',np.shape(newHERApspec)
print '# of baselines: ',np.shape(pairs)


##### Old OQE Module PSPEC #####
dsv = DataSet()
dsn = DataSet()

data_dict_v = {}
data_dict_n = {}
flg_dict = {}
evenodd = [0,1]
keys = []
for key in bls:
    for eo in evenodd:
        newkey = key +(eo,)
        data_dict_v[newkey] = d[eo].get_data(key)
        data_dict_n[newkey] = np.zeros_like(d[eo].get_data(key))
        flg_dict[newkey] = np.logical_not(d[eo].get_flags(key))
        keys.append(newkey)
conj_dict = None
dsv.set_data(dsets=data_dict_v, conj=conj_dict, wgts=flg_dict)
dsn.set_data(dsets=data_dict_n, conj=conj_dict, wgts=flg_dict)
pIv, pIn  = make_PS(keys, dsv, dsn, grouping=False)
print np.shape(pIv)
def fold(pspec):
    return np.append(pspec[10],np.mean((pspec[11:],pspec[:10][::-1]),0))

### Scale power spectra using code from Elderberry ###
freqs = np.linspace(.1,.2,203)
afreqs = freqs[95:116]
WINDOW='none'
fq = np.mean(afreqs)
z = C.pspec.f2z(fq)
sdf = freqs[1] - freqs[0]
PAPER_BEAM_POLY = [ -1.55740671e+09,  1.14162351e+09, -2.80887022e+08,  9.86929340e+06, 7.80672834e+06, -1.55085596e+06,  1.20087809e+05, -3.47520109e+03]
bm = np.polyval(PAPER_BEAM_POLY, fq) * 2.35
B = sdf * afreqs.size * C.pfb.NOISE_EQUIV_BW[WINDOW]
scalar = C.pspec.X2Y(z) * bm * B

print 'Scalar: ',scalar
print 'Shape of pIv: ',np.shape(pIv)
print 'Shape of newHERApspec: ',np.shape(newHERApspec)

# Scale the old and new PSPECs
eberryPSPEC = scalar*np.mean(np.mean(pIv,0),1)
nhoPSPEC = scalar*np.mean(np.mean(newHERApspec,0),1)

etas = np.fft.fftshift(C.pspec.f2eta(afreqs))
kpl = etas * C.pspec.dk_deta(z)
wavelength = C.cosmo_units.c/(.15*1e9)
bl_len = 30.
ubl = bl_len/wavelength
kperp = C.pspec.dk_du(.15)*ubl
k = np.sqrt(kperp**2 + kpl**2)
k[:10] *= -1.
plt.figure(figsize=(10,5))
plt.rc('text', usetex=True)
plt.semilogy(k,np.abs(nhoPSPEC.real),'b',label='New OQE')
plt.semilogy(k,np.abs(nhoPSPEC.real),'b.')
plt.semilogy(k,np.abs(eberryPSPEC.real),'r--',label='e-berry')
plt.semilogy(k,np.abs(eberryPSPEC.real),'r.')
#plt.ylabel(r'$\frac{k^3}{2 \pi}P(k) \ [mK^2]$',fontsize=30)
plt.ylabel(r'$P(k) \ [(h^{-1}Mpc)^{3} mK^2]$',fontsize=30)
plt.xlabel(r'$k \ [h \ Mpc^{-1}$]',fontsize=30)
plt.legend(loc=4,prop={'size': 20})
plt.grid()
plt.show()

