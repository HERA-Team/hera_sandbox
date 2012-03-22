"""Implements the Polyphase-Filter Bank spectral decomposition algorithm.
pfb.WINDOWS is a list of built-in windowing functions."""

# Copyright (C) 2005 Aaron Parsons
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

import numpy as n
from aipy._cephes import i0

NOISE_EQUIV_BW = {
    'blackman': 1.73,
    'blackman-harris': 2.01,
    'gaussian0.4': 1.45,
    'hamming': 1.37,
    'hanning': 1.50,
    'kaiser2': 1.50,
    'kaiser3': 1.80,
    'parzen': 1.33,
    'none': 1,
}
WINDOWS = NOISE_EQUIV_BW.keys()

pm = {}
pm['L'] = None
pm['taps'] = None
pm['fwidth'] = None
pm['window'] = None
pm['window_name'] = None
pm['sinx_x'] = None
pm['window_sinx_x'] = None

def __set_pm__(L, window, taps, fwidth):
    global pm
    if pm['L'] == L and pm['taps'] == taps and \
       pm['fwidth'] == fwidth:
        if type(window) == str and pm['window_name'] == window: return
        elif window is pm['window']: return
    else:
        pm['L'] = L
        pm['taps'] = taps
        pm['fwidth'] = fwidth
        def sinx_x(x):
            t = n.pi * taps * fwidth * (x/float(L) - .5)
            v = n.where(t != 0, t, 1)
            return n.where(t != 0, n.sin(v) / v, 1)
        pm['sinx_x'] = n.fromfunction(sinx_x, (L,))
    if type(window) == str: 
        wf = {}
        wf['blackman'] = lambda x: .42-.5*n.cos(2*n.pi*x/(L-1))+.08*n.cos(4*n.pi*x/(L-1))
        wf['blackman-harris'] = lambda x: .35875 - .48829*n.cos(2*n.pi*x/(L-1)) + .14128*n.cos(4*n.pi*x/(L-1)) - .01168*n.cos(6*n.pi*x/(L-1))
        wf['gaussian0.4'] = lambda x: n.exp(-0.5 * ((x - (L-1)/2)/(0.4 * (L-1)/2))**2)
        wf['kaiser2'] = lambda x: i0(n.pi * 2 * n.sqrt(1-(2*x/(L-1) - 1)**2)) / i0(n.pi * 2)
        wf['kaiser3'] = lambda x: i0(n.pi * 3 * n.sqrt(1-(2*x/(L-1) - 1)**2)) / i0(n.pi * 3)
        wf['hamming'] = lambda x: .54 - .46 * n.cos(2*n.pi*x/(L-1))
        wf['hanning'] = lambda x: .5 - .5 * n.cos(2*n.pi*x/(L-1))
        wf['parzen'] = lambda x: 1 - n.abs(L/2. - x) / (L/2.)
        wf['none'] = lambda x: 1
        pm['window'] = n.fromfunction(wf[window], (L,))
        pm['window_name'] = window
    else: 
        pm['window'] = window
        pm['window_name'] = None
    pm['window_sinx_x'] = pm['window'] * pm['sinx_x']

def __pfb_fir__(data, window='hamming', taps=8, fwidth=1):
    L = data.shape[-1]
    __set_pm__(L, window, taps, fwidth)
    d = data * pm['window_sinx_x']
    try: d.shape = d.shape[:-1] + (taps, L/taps)
    except: raise ValueError("More taps than samples")
    return n.sum(d, len(d.shape) - 2)

def pfb(data, window='hamming', taps=8, fwidth=1, fft=n.fft.fft):
    """Perform PFB on last dimension of 'data' for multi-dimensional arrays.
    'window' may be a name (e.g. 'hamming') or an array with length of the
    last dimension of 'data'.  'taps' is the number of PFB taps to use.  The
    number of channels out of the PFB will be length out of the last 
    dimension divided by the number of taps. 'fwidth' scales the width of 
    each channel bandpass (to create overlapping filters, for example)."""
    return fft(__pfb_fir__(data, window, taps, fwidth))

def streaming_pfb(data, nfreq, window='hamming', taps=8, fwidth=1,
        fft=n.fft.fft):
    """Perform as many PFBs of length 'nfreq' as fit in the 1d array 'data'.
    PFBs are computed using windows of length 'taps' * 'nfreq', which step
    through 'data' in increments of 'nfreq'."""
    d = n.resize(data, (n.floor(len(data)/nfreq), nfreq))
    for t in range(taps-1):
        d = n.concatenate([d[:-1], d[1:,t*nfreq:(t+1)*nfreq]], axis=1)
    return pfb(d, window=window, taps=taps, fwidth=fwidth, fft=fft)

def streaming_fft(data, nfreq, fft=n.fft.fft):
    """Perform as many FFTs of length 'nfreq' as fit in the 1d array 'data'."""
    d = n.resize(data, (n.floor(len(data)/nfreq), nfreq))
    return fft(d)
