#! /usr/bin/env python
import capo.hex as hx, capo.arp as arp, capo.red as red, capo.omni as omni
import numpy as n, aipy as a
import sys,optparse
import numpy as np

o = optparse.OptionParser()
a.scripting.add_standard_options(o,cal=True,pol=True)
o.add_option('--ubls', default='', help='Unique baselines to use, separated by commas (ex: 1_4,64_49).')
o.add_option('--ex_ants', default='', help='Antennas to exclude, separated by commas (ex: 1,4,64,49).')
o.add_option('--outpath', default=None,help='Output path of solution npz files. Default will be the same directory as the data files.')
o.add_option('--plot', action='store_true', default=False, help='Turn on plotting in firstcal class.')
o.add_option('--verbose', action='store_true', default=False, help='Turn on verbose.')
opts,args = o.parse_args(sys.argv[1:])

def flatten_reds(reds):
    freds = []
    for r in reds:
        freds += r
    return freds
#get frequencies
uv = a.miriad.UV(args[0])
fqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

#hera info assuming a hex of 19 and 128 antennas
aa = a.cal.get_aa(opts.cal, fqs)
ex_ants = []
ubls = []
for a in opts.ex_ants.split(','):
    try: ex_ants.append(int(a))
    except: pass
for bl in opts.ubls.split(','):
    try:
        i,j = bl.split('_')
        ubls.append((int(i),int(j)))
    except: pass
print 'Excluding Antennas:',ex_ants
if len(ubls) != None: print 'Using Unique Baselines:',ubls
info = omni.aa_to_info(aa, fcal=True, ubls=ubls, ex_ants=ex_ants)
reds = flatten_reds(info.get_reds())

print 'Number of redundant baselines:',len(reds)
#Read in data here.
ant_string =','.join(map(str,info.subsetant))
bl_string = ','.join(['_'.join(map(str,k)) for k in reds])

times, data, flags = arp.get_dict_of_uv_data(args, bl_string, opts.pol, verbose=True)
datapack,wgtpack = {},{}
for (i,j) in data.keys():
    datapack[(i,j)] = data[(i,j)][opts.pol]
    wgtpack[(i,j)] = np.logical_not(flags[(i,j)][opts.pol])
dlys = n.fft.fftshift(n.fft.fftfreq(fqs.size, np.diff(fqs)[0]))

#gets phase solutions per frequency.
fc = omni.FirstCal(datapack,wgtpack,fqs,info)
sols = fc.run(finetune=True,verbose=opts.verbose,plot=opts.plot,noclean=False,offset=False,average=True,window='none')

#Save solutions
if len(args)==1: filename=args[0]
else: filename='fcgains.%s.npz'%opts.pol #if averaging a bunch together of files together.
if not opts.outpath is None:
    outname='%s/%s'%(opts.outpath,filename.split('/')[-1])
else:
    outname='%s'%filename
omni.save_gains_fc(sols,fqs, opts.pol[0], outname, ubls=ubls, ex_ants=ex_ants)
