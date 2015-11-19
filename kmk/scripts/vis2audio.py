#! /usr/bin/env python
import aipy as a
import capo
import numpy as np
import matplotlib.pylab as pl
import glob
import cmath
import scipy.io.wavfile as wav

data_dir = '/home/kara/capo/kmk/data/'
data_files = sorted(glob.glob(data_dir+'*.xx.uvcRRE'))
audio = []
for file in data_files:
    uv = a.miriad.UV(file)
    N = uv.__getitem__('nchan')

    uv.select('antennae',90,91,include=True)
    time = []
    for (uvw, t, (i,j)), data, flag in uv.all(raw=True):
        bl = a.miriad.ij2bl(i,j)
        time.append(t)
        # get all data for one baseline and one frequency
        audio.append(data[N/2-1])

audio = map(cmath.phase, audio)
dt = 0.3
ref_freq = 300
ref_signal = np.arange(int(dt * 44100))
ref_signal = np.sin(2*np.pi* ref_freq / 44100 * ref_signal)
wavfile = np.array([])
for note in audio:
    wavfile = np.r_[wavfile, (note + np.pi) * ref_signal * 32750 / (2*np.pi)]
wavfile = np.array([wavfile, wavfile], dtype=np.int16).transpose()
wav.write('maybe-fringe.wav', 44100, wavfile)

pl.plot(audio)
pl.show()

