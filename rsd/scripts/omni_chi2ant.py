#!/usr/bin/env python2

################################################################################
## TODO write a sentence here.
## Copyright (C) 2014  Rachel Domagalski: domagalski@berkeley.edu
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## ## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

import argparse
import numpy as np
import omnical.info as oi
import omnical.calib as oc
import matplotlib.pyplot as plt

def get_blperant(old_info):
    info = oi.RedundantInfoLegacy() #reading old firstcal files
    info.fromfile(old_info)
    reds = info.get_reds()
    antpos = info.get_antpos()
    new_info = oi.RedundantInfo()
    new_info.init_from_reds(reds, antpos)
    blperant = new_info.blperant
    return (info.subsetant, blperant)

def get_nfreq(filenames):
    """
    This function takes in a list of files and returns the number of
    frequency channels. All files must have same number of channels.
    """
    def __nf(fname):
        omnidata = np.load(fname)
        omnikeys = filter(lambda s: ',chisq' in s, omnidata.iterkeys())
        _, nfrequencies = omnidata[omnikeys[0]].shape
        nantennas = len(omnikeys) - 1 # One of these is the total chisq
        omnidata.close()
        return (nfrequencies, nantennas)

    nfreq, nant = map(set, zip(*map(__nf, filenames)))
    if len(nfreq) != 1:
        raise ValueError('Frequency channels do not match across all data.')
    if len(nant) != 1:
        raise ValueError('Number of antennas does not match across all data.')
    return (nfreq.pop(), nant.pop())

def grab_freq(filename, subsetant, freq, pol):
    """
    Select data from a file along a certain frequency channel.
    """
    omnidata = np.load(filename)
    antchisq = np.array([omnidata[pol + ',chisq%d'%a] for a in subsetant])
    chisq = omnidata[pol + ',chisq']
    if freq == None:
        fstart = 0
        fend = chisq.shape[1]
    elif '_' in freq:
        fstart, fend = map(int, freq.split('_'))
        fend += 1
    else:
        fstart = int(freq)
        fend = fstart + 1
    freq_select = np.sum(antchisq[:,:,fstart:fend], axis=2)
    active_channels = np.sum(chisq[:,fstart:fend] != 0, axis=1)
    omnidata.close()
    return np.transpose(freq_select / active_channels)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('omni_output', nargs='+', help='Omnical outputs to be plotted.')
    parser.add_argument('-f', '--freq',
                        help='Select a frequency channel or a frequency range.')
    parser.add_argument('-i', '--info', required=True,
                        help='Redundant info file for the data being plotted.')
    parser.add_argument('-l', '--log', action='store_true',
                        help='Plot on a log scale.')
    parser.add_argument('-p', '--pol', metavar='POL',
                        default='xx', choices=['xx', 'xy', 'yy', 'XX', 'XY', 'YY'],
                        help='Polarization to plot (xx, xy, yy).')

    # Parse args and organize input info
    args = parser.parse_args()
    pol = args.pol
    freq = args.freq
    redinfo = args.info
    infiles = args.omni_output
    plot_log = args.log
    nfreq, nant = get_nfreq(infiles)

    # Select the frequency data
    chi2ant = np.empty((0, nant))
    chisq = np.array([])
    subsetant, blperant = get_blperant(redinfo)

    # Get the chisq/ant data and throw it into an array. I'm also collecting the
    # chisq data to verify that the sum of the antenna chi^2 is twice the chisq
    # that gets computed, as it should be. This is in a unit test, but since
    # there's been some code that's been buggy, I'm doing it here again.
    for fname in infiles:
        chi2ant = np.append(chi2ant, grab_freq(fname, subsetant, freq, pol), axis=0)

    # Sum everything up
    sum_chi2ant = np.sum(chi2ant, axis=1)

    # Make plots.
    if plot_log:
        plt.imshow(np.log10(chi2ant/blperant), interpolation='nearest')
    else:
        plt.imshow(chi2ant/blperant, interpolation='nearest')
    plt.xlabel('Antenna index')
    plt.ylabel('Integration')
    plt.title('$\chi^2$ for each antenna')
    plt.colorbar()
    plt.show()
