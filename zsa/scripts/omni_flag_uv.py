#! /usr/bin/env python
import numpy as n
import aipy as a, capo as C
import sys,glob,os
#import omnical.calibration_omi as omni

#Run on day to day basis.

#yea....dont actually need this 
def getbl2omni():
    bl2omni = {}
    info = omni.read_redundantinfo('redundantinfo_PSA64_ba19_37_50.bin')
    antsubset = info['subsetant']
    ubllist = info['ublindex']
    for sep in xrange(len(ubllist)):
        for num in xrange(len(ubllist[sep])):
            i,j = ubllist[sep][num][0:2]
            ai = antsubset[int(i)]
            aj = antsubset[(j)]
            bl2omni[(ai,aj)] = sep
    return bl2omni, len(ubllist)
            
#bl2omni, nseps = getbl2omni()
args = sys.argv[1:]

#get sigmas. Since omnical files are chunks. Read all uv files in that chunk 
#all we want the omnical files for is the chi_2. Get that first.
uvtimes, uvdata, uvflags = C.arp.get_dict_of_uv_data(args, antstr='cross', polstr='I', verbose=True)
npzfiles = {}


jday = args[0].split('.')[1]
print jday
sig = {}
pols = ['I']
for p in ['xx','yy']:
    npzfiles[p] = glob.glob('/home/aparsons/omnical_results/pass02/data_psa*_%s.*%s*v2.npz'%(jday,p))

#want to get chi2 for stokes chi2I = .5(chi2xx + chi2yy)
#print npzfiles

ichi2 = {}
for p in ['xx','yy']:
    #get chi_2 from omnical files. These are for the entire array time vs fq.
    #get it for a whole chunk XXX make npz files
    ichi2[p] = []
    for npzs in npzfiles[p]:
        npz = n.load(npzs)
        ichi2[p].append(n.where(n.isnan(npz['chi2_lin']), 0, 1./npz['chi2_lin']))
    ichi2[p] = n.concatenate(ichi2[p]) #should be time ordered. This 

print 'Combining pols to make stokes I chi2'
ichi2i = .5*(ichi2['xx'] + ichi2['yy'])

for p in pols:
    sig[p] = {}
    for bl in uvdata.keys():
        d,w = [],[] 
        w = n.logical_not(uvflags[bl][p])*ichi2i #flags is True when data flagged.
        d = uvdata[bl][p]
        w = n.where(n.isnan(d), 0, w)
        w = n.where(w< 1e-8, 0, w)
        d = n.where(w>0, d, 0)
        sig[p][bl] = n.sqrt(n.median(n.abs(d)**2, axis=0))
        
print sig[p].keys()

def mfunc(uv, p, d, f):
    uvw,t,(i,j) = p
    if i==j:
        return p,d,f
    bl = a.miriad.ij2bl(i,j)
    pol = a.miriad.pol2str[uv['pol']]
    f = n.logical_or( n.where(n.abs(d) > 3*sig[pol][bl], True, False), f) #Find where data is above 3 sigma for that baseline. OR with incoming flags.
    d = n.where(f, 0, d)# where flags==True make 0, otherwise, data.
    return p, d, f


for filename in args:
    outfile=os.getcwd()+'/'+filename.split('/')[-1]+'f'
    if os.path.exists(outfile):
        print '    %s exists.  Skipping...' % outfile
        continue
    print '     Writing', outfile
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, append2hist=' '.join(sys.argv)+'\n', raw=True)
    

