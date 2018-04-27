"""
fourier_interp.py
-----------------

Fourier Interpolation
"""
import numpy as np
import os
import sys
import copy
from scipy import signal
from collections import OrderedDict as odict
from hera_cal.datacontainer import DataContainer


def fourier_filter(vis, flags, kernel_width=10, kernel='tophat', axis=1, stop_tol=1e-2, maxiter=5, 
                  copy_vis=True, edge_cut=0):
    """
    Fast fourier filtering + interpolation for flagged data.

    Parameters:
    -----------
    vis : type=complex ndarray, holding visibility data, with [0] indexing time and [1] freq

    flags : type=boolean ndarray, holding visibility flags, matching shape of vis

    kernel_width : type=int, width of kernel in bin number. If axis=1, this represents
                   the number of delay modes, if axis=0 this represents number of fringe-rate modes.

    kernel : type=str, type of kernel for smoothing, options=['tophat', 'gaussian']

    axis : type=int, axis along vis to perform smoothing

    stop_tol : type=float, stopping tolerance for median of residuals

    maxiter : type=int, max number of iterations

    copy_vis : type=boolean, if True, make a copy of vis so as not to overwrite input reference

    edge_cut : type=int, number of channels to ignore on either side of band edge

    Output: (vis, vis_hat)
    -------
    vis : type=complex ndarray, holding visibility data with RFI flags filled-in
    
    vis_hat : type=complex ndarray, holding final smoothed estimate of the visibility
    """
    # copy
    if copy_vis:
        vis = copy.copy(vis)
        flags = copy.copy(flags)

    # get vis shape
    Ntimes = vis.shape[0]
    Nfreqs = vis.shape[1]

    # get FFT kernel
    if kernel == 'tophat':
        # setup tophat kernel
        kern = signal.windows.boxcar(vis.shape[axis])
        if axis == 0:
            kern[kernel_width:Ntimes-kernel_width] = 0.0
            kern = np.repeat(kern[:, None], Nfreqs, axis=1)
        elif axis == 1:
            kern[kernel_width:Nfreqs-kernel_width] = 0.0
            kern = np.repeat(kern[None, :], Ntimes, axis=0)

    elif kernel == 'gaussian':
        # setup Gaussian kernel
        kern = np.fft.fftshift(signal.windows.gaussian(vis.shape[axis], kernel_width))
        if axis == 0:
            kern = np.repeat(kern[:, None], Nfreqs, axis=1)
        elif axis ==1:
            kern = np.repeat(kern[None, :], Ntimes, axis=0)

    else:
        raise ValueError("didn't recognize {} kernel".format(kernel))

    # null flagged pixels
    vis[flags] *= 0

    # null band edges
    vis[:, :edge_cut] *= 0.0
    vis[:, -edge_cut:] *= 0.0
    flags[:, :edge_cut] = False
    flags[:, -edge_cut:] = False

    # iterate over stopping tolerance
    med_residual = 100.0
    niters = 0
    while med_residual > stop_tol:
        # FFT w/ window along axis
        vis_fft = np.fft.fft(vis, axis=axis)

        # multiply by kernel
        vis_fft *= kern

        # transform back
        vis_hat = np.fft.ifft(vis_fft, axis=axis)

        # estimate residuals
        med_residual = np.median(np.abs(np.abs(vis[flags]) - np.abs(vis_hat[flags])))

        # fill in RFI flags
        vis[flags] = vis_hat[flags]

        # break reached stopping tol
        if med_residual <= stop_tol or niters >= maxiter:
            break

        niters += 1

    return vis, vis_hat


def fourier_interp(data, flags, kernel_width=10, kernel='tophat', axis=1, stop_tol=1e-2, maxiter=5, 
                   copy_vis=True):
    """
    Fourier filtering + interpolation for flagged data

    Parameters:
    -----------
    data : type=DataContainer, holding complex visibility data

    flags : type=DataContainer, holding boolean flag arrays for data

    Keyword parameters passed to fourier_filter()

    Output: (interp_data, model)
    -------
    interp_data : type=DataContainer, holding complex visibility data with
                  flags filled in with model

    model : type=DataContainer, holding smoothed visibility model
    """
    # setup output data structures
    interp_data = odict()
    model = odict()

    # iterate over keys
    for k in data.keys():
        # perform fourier fit
        d_int, mdl = fourier_filter(data[k], flags[k], kernel_width=kernel_width, kernel=kernel,
                                    axis=axis, stop_tol=stop_tol, maxiter=maxiter, copy_vis=copy_vis)

        # insert into structures
        interp_data[k] = d_int
        model[k] = mdl

    interp_data = DataContainer(interp_data)
    model = DataContainer(model)

    return interp_data, model


