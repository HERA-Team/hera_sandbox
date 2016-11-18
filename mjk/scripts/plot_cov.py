#! /usr/bin/env python
import matplotlib
#matplotlib.use('Agg')
import aipy as a, numpy as n, pylab as p, capo
import glob, optparse, sys, random
from IPython import embed
o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, cal=True)
o.add_option('--plot', action='store_true',
    help='Generate plots')
o.add_option('--window', dest='window', default='blackman-harris',
    help='Windowing function to use in delay transform.  Default is blackman-harris.  Options are: ' + ', '.join(a.dsp.WINDOW_FUNC.keys()))
o.add_option('--output', type='string', default='',
    help='output directory for pspec_boot files (default "")')
o.add_option('--cav',action='store_true',
    help='Use Cav to inplace of Covariance Matrix')
o.add_option('--auto', action='store_true',
    help='Scales main diagonal of Cav by auto-covariances')

opts,args = o.parse_args(sys.argv[1:])

random.seed(0)
POL = 'I'
LST_STATS = False
DELAY = False
NGPS = 5
INJECT_SIG = 0.
SAMPLE_WITH_REPLACEMENT = True
NOISE = .0
PLOT = opts.plot
def get_data(filenames, antstr, polstr, rmbls, verbose=False):
    # XXX could have this only pull channels of interest to save memory
    lsts, dat, flg = [], {}, {}
    if type(filenames) == 'str': filenames = [filenames]
    for filename in filenames:
        if verbose: print '   Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, antstr, polstr)
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            bl = a.miriad.ij2bl(i,j)
            if bl in rmbls: continue
            lst = uv['lst']
            if len(lsts) == 0 or lst != lsts[-1]: lsts.append(lst)
            if not dat.has_key(bl): dat[bl],flg[bl] = [],[]
            dat[bl].append(d)
            flg[bl].append(f)
            #if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
            #pol = a.miriad.pol2str[uv['pol']]
            #if not dat[bl].has_key(pol):
            #    dat[bl][pol],flg[bl][pol] = [],[]
            #dat[bl][pol].append(d)
            #flg[bl][pol].append(f)
    return n.array(lsts), dat, flg

def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X = n.array(m, ndmin=2, dtype=n.complex)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1)
    return (n.dot(X, X.T.conj()) / fact).squeeze()

def noise(size):
    return n.random.normal(size=size) * n.exp(1j*n.random.uniform(0,2*n.pi,size=size))


def get_cav(c_mat,nchan,scaling=False):
        #Compute Cav by taking diags, averaging and reforming the matrix
        diags =[c_mat.diagonal(count) for count in xrange(nchan-1, -nchan,-1)]
        cav=n.zeros_like(c_mat)

        for count,count_chan in enumerate(range(nchan-1,-nchan,-1)):
                cav += n.diagflat( n.mean(diags[count]).repeat(len(diags[count])), count_chan)

        if scaling:
            temp=cav.copy()
            for count in xrange(nchan):
                cav[count,:] *= n.sqrt(c_mat[count,count]/temp[count,count])
                cav[:,count] *= n.sqrt(c_mat[count,count]/temp[count,count])

        if not n.allclose(cav.T.conj(),cav):
            print 'Cav is not Hermitian'
            print 'bl'

        if PLOT and False:
            p.subplot(131); capo.arp.waterfall(c_mat,mode='real',mx=3,drng=6);
            p.title('C')
            p.subplot(132); capo.arp.waterfall(cav,mode='real',mx=3,drng=6); p.colorbar()
            p.title('Cav')
            p.subplot(133); capo.arp.waterfall(diff,mode='real',mx=500,drng=500); p.colorbar()
            p.title('Abs Difference')
            p.suptitle('%d_%d'%a.miriad.bl2ij(bl))
            p.show()

        return cav

def write_cav_stats(data,nchan):
    if opts.auto:
        statfile='Cav_{0}_data_{1}_scaled.csv'.format(k,opts.sep)
        picfile='All_cav_diff_{0}_{1}_scaled.png'.format(k,opts.sep)
    else:
        statfile='Cav_{0}_data_{1}.csv'.format(k)
        picfile='All_cav_diff_{0}_{1}.png'.format(k)
    if not opts.output == '':
        statfile =opts.output+'/'+statfile
        picfile =opts.output+'/'+picfile
    f = open(statfile,'w')
    f.write('Baseline \t Max \t Min \t Mean \t Median \t Eig \t Eig_Cav \t Rms_C_Cav_diff\n')
    fig =p.figure()

    for num,bl in enumerate(data):
        c_mat=cov(data[bl])
        cav = get_cav(c_mat,nchan,scaling=opts.auto)
        diff=abs(c_mat-cav)
        U_cav,S_cav,V_cav=n.linalg.svd(cav.conj())
        U_c,S_c,V_c=n.linalg.svd(c_mat.conj())


        if PLOT and False
            f=p.figure()
            p.plot(S_c,label='C')
            p.plot(S_cav,label='Cav')
            p.yscale('log')
            p.legend(loc='best')
            p.title('Eigenvalues')
            p.show()
        ###Plot and making file fore residuals, amplitudes of eigenvalues, and variance of data

        f.write(' {0} \t {1} \t {2} \t {3} \t {4} \t {5} \t {6} \t {7} \n'.format(bl, n.max(diff), n.min(diff),n.mean(diff), n.median(diff), S_c[0],S_cav[0], n.sqrt(n.mean(diff**2))))

        p.subplot(7,8,num+1); capo.arp.waterfall(diff/n.sqrt(nchan),mode='real',mx=100,drng=100); p.title(bl); p.axis('off')

    f.close()
    fig.subplots_adjust(hspace=.6)
    p.suptitle('Percent diff/sqrt({0})'.format(nchan))
    fig.savefig(picfile)
    p.close()

def plot_eig(data,nchan):
    days=data.keys()
    eig_order=[]
    eigs = []
    eigs_cav = []
    for bl in data:
       c_mat=cov(data[bl])
       cav = get_cav(c_mat,nchan,scaling=opts.auto)
       U,S,V= n.linalg.svd(c_mat.conj())
       U_cav,S_cav,V_cav = n.linalg.svd(cav.conj())
       eig_order.append(S[0])
       eigs.append( n.fft.fftshift(n.fft.fft(V.T.conj(),axis=0)))
       eigs_cav.append( n.fft.fftshift(n.fft.fft(V_cav.T.conj(),axis=0)))


    order=n.argsort(eig_order)

    eig_order=n.take(eig_order,order)
    eigs=n.take(eigs,order,axis=0)
    eigs_cav=n.take(eigs_cav,order,axis=0)
    embed()
    fig=p.figure(1)
    for cnt,eig in enumerate(eigs):
        p.plot(eig[0] + cnt*5)
    p.title('Eigenvectors for day {0}'.format(k))
    p.show()
    p.savefig('eigenvectors_{0}.png'.format(k))
    p.clf()
    for cnt,eig in enumerate(eigs_cav):
        p.plot(eig[0] + cnt*5)
    p.title('Eigenvectors of Cav for day {0}'.format(k))
    p.savefig('eigenvectors_cav_{0}.png'.format(k))
    p.clf()
    p.close()


dsets = args


WINDOW = opts.window
uv = a.miriad.UV(dsets.values()[0][0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
del(uv)

afreqs = freqs.take(chans)
nchan = chans.size

aa = a.cal.get_aa(opts.cal, n.array([.150]))
bls,conj = capo.red.group_redundant_bls(aa.ant_layout)
jy2T = capo.pspec.jy2T(afreqs)
window = a.dsp.gen_window(nchan, WINDOW)
#if not WINDOW == 'none': window.shape=(1,nchan)
if not WINDOW == 'none': window.shape=(nchan,1)

# acquire the data
#antstr = '41_49,3_10,9_58,22_61,20_63,2_43,21_53,31_45,41_47,3_25,1_58,35_61,42_63,2_33'
antstr = 'cross'
lsts,data,flgs = get_data(dsets, antstr=antstr, polstr=POL, rmbls=rmbls, verbose=True)

x={}
for bl in data:
    print k, bl
    d = data[bl][:,chans] * jy2T
    if conj[bl]: d = n.conj(d)
    x[bl] = n.transpose(d, [1,0]) # swap time and freq axes
        
bls_master = x.values()[0].keys()
nbls = len(bls_master)
print 'Baselines:', nbls

##write file of stats for c and cav
write_cav_stats(x,nchan)
plot_eig(x,nchan)
# Compute baseline auto-covariances and apply inverse to data
_Ix = {}
C,_C,_Cx = {},{},{}
Cav = {}
for bl in x:
    C[bl] = cov(x[bl])

    if opts.cav:
        C[bl]=get_cav(C[bl],nchan,scaling=opts.auto)

    I[bl] = n.identity(C[bl].shape[0])
    U,S,V = n.linalg.svd(C[bl].conj())

    if S[0] >= 18:
        C[bl]=cov(x[bl])
        U,S,V=n.linalg.svd(C[bl].conj())

    _C[bl] = n.einsum('ij,j,jk', V.T, 1./S, U.T)
    _Cx[bl] = n.dot(_C[bl], x[bl])
    _Ix[bl] = x[bl].copy()


    if PLOT and True:
        #p.plot(S); p.show()
        p.subplot(311); capo.arp.waterfall(x[bl], mode='real')
        p.subplot(334); capo.arp.waterfall(C[bl])
        p.subplot(335); p.plot(n.einsum('ij,jk',n.diag(S),V).T.real)
        p.subplot(336); capo.arp.waterfall(_C[bl])
        p.subplot(313); capo.arp.waterfall(_Cx[bl], mode='real')
        p.suptitle('%d_%d'%a.miriad.bl2ij(bl))
         p.figure(2); p.plot(n.diag(S))
        p.show()

