#!/usr/bin/env python2.7
"""
smooth_calfits.py
---------------

load calfits gain solutions,
and then smooth them using median filtering,
polynomial fitting or GP fitting.

Nicholas Kern
Dec. 2017
"""
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from pyuvdata import UVCal, UVData
import numpy as np
import argparse
import os
import scipy.signal as signal
from sklearn import gaussian_process as gp
import copy
from make_calfits import make_calfits
from sklearn import linear_model
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures

args = argparse.ArgumentParser(description="")

# Required Parameters
args.add_argument("files", type=str, nargs='*', help='path to calfits file(s) to smooth')
args.add_argument("--ext", type=str, default='smooth', help="output calfits file extension")
# IO Parameters
args.add_argument("--outdir", default=None, type=str, help="output directory")
args.add_argument("--overwrite", default=False, action='store_true', help='overwrite output files')
args.add_argument("--silence", default=False, action='store_true', help='silence output to stdout')
# Medfilt Parameters
args.add_argument("--medfilt", default=False, action='store_true', help="perform median filtering across freq on real and imag components")
args.add_argument("--medfilt_kernel", default=21, type=int, help='size of median filter kernel in freq channels')
args.add_argument("--medfilt_flag", default=False, action='store_true', help='use median filtering to flag channels with large amplitude zscores')
# Polyfit Parameters
args.add_argument("--polyfit", default=False, action='store_true', help='fit a polynomial across frequency')
args.add_argument("--polyRANSAC", default=False, action='store_true', help='instead of linear regressor, use RANSAC regressor')
args.add_argument("--polydeg", default=2, type=int, help='polynomial degree across time and frequency')
# Gaussian Process Parameters
args.add_argument("--gpfit", default=False, action='store_true', help='fit a gaussian process to real and imag comps of data. recommended to precede with median filter.')
args.add_argument("--gp_max_dly", default=1e2, type=float, help="maximum delay in nanosecs of GP length scale hyper parameter")
args.add_argument("--gp_dly", default=1e2, type=float, help="starting point for delay hyperparameter")
args.add_argument("--gp_min_min", default=2, type=float, help="minimum length scale across time in minutes")
args.add_argument("--gp_min", default=10, type=float, help="starting point for minute hyperparameter")
args.add_argument("--gp_nrestart", default=1, type=int, help="number of hyperparameter restarts")
args.add_argument("--gp_freq_thin", default=8, type=int, help="thinning factor across frequency before GP fit to data")
args.add_argument("--gp_time_thin", default=2, type=int, help="thinning factor across time before GP fit to data")
args.add_argument("--gp_meanfunc", default=False, action='store_true', help="fit a polynomial to use as mean function before GP fit, uses poly_freqdeg and poly_timedeg arguments")
args.add_argument("--gp_avgtime", default=False, action='store_true', help='average across time before fitting GP')
args.add_argument("--gp_optimizer", default='fmin_l_bfgs_b', help="GP optimizer")
def echo(message, type=0, verbose=True):
    if verbose:
        if type == 0:
            print(message)
        elif type == 1:
            print('\n{}\n{}'.format(message, '-'*40))

def smooth_calfits(file, ext='smooth', outdir=None, overwrite=False,
                  medfilt=False, medfilt_kernel=7, medfilt_flag=False,
                  polyfit=False, polydeg=2, polyRANSAC=False,
                  gpfit=False, gp_max_dly=1e2, gp_nrestart=1, gp_freq_thin=8, gp_time_thin=2,
                  gp_meanfunc=False, gp_avgtime=False, gp_flagrm=False, gp_min_min=5.0,
                  gp_dly=1e2, gp_min=10, gp_optimizer='fmin_l_bfgs_b', verbose=True):
    """
    """
    # construct output file
    if outdir is None:
        outdir = os.path.dirname(file)

    output_fname = os.path.splitext(os.path.basename(file))
    output_fname = os.path.join(outdir, output_fname[0] + ".{}".format(ext) + output_fname[1])

    # check overwrite
    if os.path.exists(output_fname) and overwrite is False:
        raise IOError("{} exists, not overwriting".format(output_fname))

    # load gains
    uvc = UVCal()
    uvc.read_calfits(file)

    # check gain array
    if uvc.cal_type != "gain":
        raise AttributeError("only works on cal_type == 'gain'")

    # extract gains and flags
    gains = uvc.gain_array
    flags = uvc.flag_array

    # get data params
    Ntimes = uvc.Ntimes
    Nfreqs = uvc.Nfreqs
    Njones = uvc.Njones
    Nspws = uvc.Nspws
    ants = uvc.ant_array
    pols = uvc.jones_array
    spws = uvc.spw_array
    freqs = uvc.freq_array.squeeze()
    time_array = uvc.time_array * 24 * 60
    time_array -= time_array.min()

    if gp_optimizer == "None":
        gp_optimizer = None

    # iterate over antennas
    for ii, a in enumerate(ants):
        echo("...working on antenna {}".format(a), type=1, verbose=verbose)

        # skip completely flagged antennas
        if np.min(uvc.flag_array[ii]):
            continue

        # iterate over polarizations
        for jj, p in enumerate(pols):

            # iterate over spws
            for kk, s in enumerate(spws):

                # get antenna gains
                ant_gains = copy.copy(gains[ii, kk, :, :, jj])
                ant_flags = copy.copy(flags[ii, kk, :, :, jj])

                ###################
                ## Median Filter ##
                ###################
                if medfilt or medfilt_flag:
                    echo("...median filtering", verbose=verbose)
                    real_medfilt = signal.medfilt(ant_gains.real, kernel_size=(medfilt_kernel, 1))
                    imag_medfilt = 1j * signal.medfilt(ant_gains.imag, kernel_size=(medfilt_kernel, 1))
                    medfilt_gains = real_medfilt + imag_medfilt

                    # flag data
                    if medfilt_flag:
                        # get residual and MAD from unfiltered and filtered data
                        residual = (np.abs(ant_gains) - np.abs(medfilt_gains))
                        residual[ant_flags] *= np.nan
                        residual = residual.squeeze()
                        resid_std = np.nanmedian(np.abs(residual - np.nanmedian(residual, axis=1)[:, np.newaxis]), axis=1) * 1.5

                        # identify outliers as greater than 10-sigma
                        bad = np.array([np.abs(residual[i]) > resid_std[i]*10 for i in range(residual.shape[0])])

                        # add to flags
                        ant_flags += bad

                    # save medfilt gains
                    if medfilt:
                        ant_gains = medfilt_gains

                ####################
                ## Polynomial Fit ##
                ####################
                if polyfit or (gp_meanfunc is True and gpfit is True):
                    echo("...fitting polynomial", verbose=verbose)
                    # setup X axis data
                    data_shape = (uvc.Nfreqs, uvc.Ntimes)
                    X = np.meshgrid(np.arange(uvc.Ntimes), np.arange(uvc.Nfreqs))
                    X_fit = np.array([X[0].ravel(), X[1].ravel()]).T.astype(np.float)
                    X_fit_thin = X_fit[~ant_flags.ravel(), :]

                    # setup y axis data
                    y_real = ant_gains.real.ravel()
                    y_real_thin = y_real[~ant_flags.ravel()]
                    y_imag = ant_gains.imag.ravel()
                    y_imag_thin = y_imag[~ant_flags.ravel()]

                    # setup polynomial model
                    if polyRANSAC:
                        regressor = linear_model.RANSACRegressor()
                    else:
                        regressor = linear_model.LinearRegression()
                    model_real = make_pipeline(PolynomialFeatures(degree=polydeg),
                                               regressor)
                    model_imag = make_pipeline(PolynomialFeatures(degree=polydeg),
                                               regressor)

                    # fit models to thinned data
                    model_real.fit(X_fit_thin, y_real_thin)
                    model_imag.fit(X_fit_thin, y_imag_thin)

                    # predict models across all data
                    y_real_pred = model_real.predict(X_fit)
                    y_imag_pred = model_imag.predict(X_fit)
                    y_pred = y_real_pred + 1j * y_imag_pred
                    y_pred = y_pred.reshape(data_shape)

                    if gpfit is False:
                        ant_gains = y_pred

                ##########################
                ## Gaussian Process Fit ##
                ##########################
                if gpfit:
                    echo("...fitting gaussian process", verbose=verbose)

                    # subtract mean function if desired
                    if gp_meanfunc:
                        ant_gains -= y_pred

                    # configure X data
                    gp_X = np.meshgrid(time_array, freqs / 1e6)
                    gp_X = np.array([gp_X[0], gp_X[1]])
                    data_shape = gp_X.shape

                    # configure y data
                    gp_yreal = ant_gains.real
                    gp_yimag = ant_gains.imag
                    gp_ydata = np.array([gp_yreal, gp_yimag])

                    # eliminate freqs when over 90% of times are flagged
                    if gp_flagrm:
                        flag_rate = np.sum(ant_flags.astype(np.float), axis=1) / Ntimes
                        good = flag_rate < 0.90
                        gp_X_thin = gp_X[:, good, :]
                        gp_ydata_thin = gp_ydata[:, good, :]
                    else:
                        gp_X_thin = gp_X
                        gp_ydata_thin = gp_ydata

                    # average across time
                    if gp_avgtime:
                        gp_X_thin = gp_X_thin[1, :, :1]
                        gp_ydata_thin = np.mean(gp_ydata_thin, axis=2)

                        # thin data
                        gp_X_thin = gp_X_thin[::gp_freq_thin, :]
                        gp_ydata_thin = gp_ydata_thin[:, ::gp_freq_thin]

                        # setup GP kernel
                        freq_lambda = 1. / (gp_max_dly * 1e-3) # MHz
                        kernel = 1**2 * gp.kernels.RBF(1./(gp_dly*1e-3), (freq_lambda, 200.0)) + gp.kernels.WhiteKernel(1e-4, (1e-8, 1e0))
                        GP = gp.GaussianProcessRegressor(kernel=kernel, optimizer=gp_optimizer, n_restarts_optimizer=gp_nrestart)

                        # get data shape
                        data_shape = data_shape[:-1]
                        gp_X = gp_X[1, :, :1]

                        # fit GP on thinned data
                        GP.fit(gp_X_thin, gp_ydata_thin.T)

                        # predict across freq array
                        gp_ypred = GP.predict(gp_X)

                        # convert to complex
                        gp_ypred = gp_ypred[:, 0] + 1j * gp_ypred[:, 1]

                        # replace gains
                        ant_gains = np.repeat(gp_ypred.reshape(-1, 1), Ntimes, axis=1)

                    else:
                        # thin data
                        gp_X_thin = gp_X_thin[:, ::gp_freq_thin, ::gp_time_thin]
                        gp_ydata_thin = gp_ydata_thin[:, ::gp_freq_thin, ::gp_time_thin]

                        # setup GP kernel
                        freq_lambda = 1. / (gp_max_dly * 1e-3) # MHz
                        kernel = 1**2 * gp.kernels.RBF((gp_min, 1/(gp_dly*1e-3)), ((gp_min_min, 20), (freq_lambda, 200.0))) + gp.kernels.WhiteKernel(1e-4, (1e-8, 1e0))
                        GP = gp.GaussianProcessRegressor(kernel=kernel, optimizer=gp_optimizer, n_restarts_optimizer=gp_nrestart)

                        # fit GP on thinned data
                        GP.fit(gp_X_thin.reshape(2, -1).T, gp_ydata_thin.reshape(2, -1).T)

                        # predict across freq array
                        gp_ypred = GP.predict(gp_X.reshape(2, -1).T).T
                        gp_ypred = gp_ypred.reshape(data_shape)

                        # convert to complex
                        gp_ypred = gp_ypred[0] + 1j * gp_ypred[1]

                        # replace gains
                        ant_gains = gp_ypred

                    # add in mean function
                    if gp_meanfunc:
                        ant_gains += y_pred

                # Assign gains
                uvc.gain_array[ii, kk, :, :, jj] = ant_gains
                uvc.flag_array[ii, kk, :, :, jj] = ant_flags

    # write to file
    echo("...saving {}".format(output_fname), type=1, verbose=verbose)
    uvc.write_calfits(output_fname, clobber=True)

if __name__ == "__main__":
    # parse args
    a = args.parse_args()

    kwargs = dict(vars(a))
    kwargs.pop("files")
    kwargs['verbose'] = kwargs['silence'] is False
    kwargs.pop('silence')
    # iterate over files
    for iii, fff in enumerate(a.files):
        smooth_calfits(fff, **kwargs)

