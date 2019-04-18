'''
This script preprocesses autocorrelations from h2c.
Given a list of observations
---------------------------------------------------
1) it runs xrfi on largest manageable chunk of data in files from a directory.
2) writes intermediary flagged miriad files for each chunk (autocorrs only)
3) concatenates all chunks into a single file.
4) delete intermediary outputs
----------------------------------------------------
'''

import numpy
from hera_qm import xrfi
from pyuvdata import UVData
import argparse
import os
cwd = os.getcwd()

parser = argparse.ArgumentParser(description='Generate rfi flagged autocorrelations using all available data at each time.')
parser.add_argument('-d','--directory',dest='directory')
#parser.add_argument('-f','--filelist',dest='filelist')
parser.add_argument('-c','--chunks',dest='chunks',default=1)
parser.add_argument('--outputdir',dest='outputdir',default=cwd)
parser.add_argument('--tempdir','--temp',dest='tempdir',default=cwd)
parser.add_argument('--output','-o',dest='output')
parser = parser.parse_args()
#user can specify a list of files or
#assert not and(parser.filelist, parser.directory)

chunk_size = parser.chunks
dir = parser.directory
tempdir = parser.tempdir
outdir = parser.outputdir
output = parser.output

file_list = glob.glob(dir+'/*.HH.uvh5')
jds = []
for fname in files:
    fname_jd = fname.split('/')[-1]
    jd = float(re.findall('[0-9]{7}.[0-9]{5}',fname_jd)[0])
    jds.append(jd)
jds = np.array(jds)
files = list(np.array(files)[np.argsort(jds)])
jds = jds[np.argsort(jds)]
files = list(np.array(files)[np.argsort(jds)])
#get diff files
files_diff = []
for jd in jds:
    files_diff.append(glob.glob(dir+'/*%s*.HH.diff.uvh5'%str(jd))[0])
chunks = int(len(files)/chunk_size)
file_chunks = [files[cnum*chunk_size:(cnum+1)*chunk_size] for cnum in range(chunks)]
file_chunks_diff = [files_diff[cnum*chunk_size:(cnum+1)*chunk_size] for cnum in range(chunks)]
jd_chunks = [jds[cnum*chunk_size:(cnum+1)*chunk_size] for cnum in range(chunks)]

chunk_files = []
chunk_files_diff = []

for chunk,chunk_d in zip(file_chunks,file_chunks_diff):
    uv = UVData()
    uv.read(chunk)
    uvd = UVData()
    uvd.read(chunk_d)
    # first round of flagging
    uvf_m, uvf_f = xrfi.xrfi_pipe(uv)
    # apply first round
    xrfi.flag_apply(uvf_f, uv, keep_existing=True, force_pol=True)
    # Second round
    uvf_m2, uvf_f2 = xrfi.xrfi_pipe(uv, alg='detrend_meanfilt')
    # Threshold
    uvf_temp = uvf_f2.copy()
    uvf_temp.to_metric(convert_wgts=True)
    uvf_final = xrfi.flag(uvf_temp, nsig_p=1.0, nsig_f=freq_threshold, nsig_t=time_threshold)
    xrfi.flag_apply(uvf_final,uv)
    xrfi.flag_apply(uvf_final,uvd)
    uv.select(bls = [(a,a) for a in np.unique(uvf_final.antenna1_array)], inplace=True)
    uvd.select(bls = [(a,a) for a in np.unique(uvf_final.antenna1_array)], inplace=True)
    chunk_files.append(tempdir+'/temp.%d.HH'%chunk)
    chunk_files_diff.append(tempdir+'/temp.diff.%d.HH'%chunk)
    uv.write_miriad([-1])
    uvd.write_miriad(chunk_files_diff[-1])

uv = UVData()
uvd = UVData()
uv.read(chunk_files)
uvd.read(chunk_files_diff)
uv.write_miriad(outputdir+output+'.HH')
uv.write_miriad(outputdir+output+'.diff.HH')
#clean up temp files
for temp,temp_d in zip(chunk_files,chunk_files_diff):
    os.rmdir(temp)
    os.rmdir(temp_d)
