import tables as tb
import numpy as np
import capo, aipy as a
import os, sys, glob

cwd = os.getcwd()
if cwd.startswith('/Users/yunfanzhang/'):
    dataDIR = '/Users/yunfanzhang/local/DATA128/DATA/'
elif cwd.startswith('/Users/yunfanz/'):
    dataDIR = '/Users/yunfanz/Data/PAPER128/DATA/'
elif cwd.startswith('/home/yunfanz/'):
    dataDIR = '/home/yunfanz/EoR/DATA128/DATA/'
sets = {
    'day1' : glob.glob(dataDIR+'zen.2456715.52*.xx.npz'),
    'day2' : glob.glob(dataDIR+'zen.2456716.52*.xx.npz'),
}
data,wgts = {}, {}
lsts = {}
chisqs = {}
for s in sets:
    if not lsts.has_key(s):
        m, g, v, x = capo.omni.from_npz(sets[s], pols='xx',verbose=True)
        lsts[s] = m['lsts']
    # #chisqs[s] = meta['chisq'][:,CH0:CH0+NCHAN]
    # for pol in vismdl:
    #     for bl in SEPS:
    #         k = (s,pol,bl)
    #         data[k] = vismdl[pol][bl][:,CH0:CH0+NCHAN]
    #         if bl in CONJ: data[k] = data[k].conj()
    #         data[k] *= bandpass[:,CH0:CH0+NCHAN]
    #         wgts[k] = np.where(np.abs(data[k]) == 0, 0., 1)
#lass Vis(IsDescription):
Vis = {}
for bls in v['xx'].keys():
	Vis[bls] = tb.ComplexAtom(itemsize=8)

import IPython; IPython.embed()
def to_hdf5(m,g,v,x,out_name='test.h5'):
	h5file = open_file(out_name, mode = "w", title = "Test file")
	meta = h5file.create_group("/", 'meta', 'Metadata including history, jds, lsts, and various chisq')
	gains = h5file.create_group("/", 'gains', 'Gains of antennas')
	vismdl = h5file.create_group("/", 'vismdl', 'Visibilities of baselines given polarization')
	xtalk = h5file.create_group("/", 'xtalk', 'cross-talk between antennas given polarization')


	for pol in ['xx']:
		table = h5file.create_table(root.vismdl, pol, Vis, "Polarization: "+pol)




