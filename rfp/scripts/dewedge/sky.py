"""Simulate sky-locked signals."""

import numpy as np
from scipy import stats


def pntsrc_foreground(
    pix,
    freqs,
    mfreq=0.15,
    Nsrcs=1000,
    intensity_bounds=(1e-3, 10),
    intensity_slope=-2,
    spectral_index_mean=-1,
    spectral_index_std=0.5,
    seed=None,
    dtype=np.float32,
):
    """
    Populate a sky with point sources; roughly ripped from hera_sim.

    Parameters
    ----------
    pix: np.ndarray
        Array with the same shape as the grid of sky pixels.
    freqs: np.ndarray of float
        Frequencies at which to evaluate the foreground signal, in GHz.
    mfreq: float, optional
        Reference frequency for adding spectral structure. Default is
        to use 150 MHz as a reference.
    Nsrcs: int, optional
        Number of sources to generate. Default is the lesser of either 1000
        or the number of sky pixels.
    intensity_bounds: length-2 array-like of float, optional
        Lower and upper bounds of the distribution from which to draw
        specific intensities, in units of Jy/beam.
    intensity_slope: float, optional
        Slope of the power-law distribution from which to draw flux densities.
    spectral_index_mean: float, optional
        Mean of the normal distribution from which spectral indices are drawn.
    spectral_index_std: float, optional
        Standard deviation of the normal distribution from which spectral indices
        are drawn.
    dtype: type, optional
        Data type to use for simulating the foreground intensities

    Returns
    -------
    intensities: np.ndarray of float
        Shape ``(pix.size, freqs.size)`` array of intensities.
    """
    if seed is not None:
        np.random.seed(seed)

    pix = np.atleast_1d(pix).astype(dtype)
    freqs = np.atleast_1d(freqs).astype(dtype)
    Nsrcs = min(Nsrcs, pix.size)
    Imin, Imax = intensity_bounds
    powerlaw_slope = -(intensity_slope + 1)
    # XXX This isn't super stable for slopes with abs < 1, but should be fine...
    intensity_scales = stats.powerlaw.rvs(
        powerlaw_slope, 1/Imax, 1/Imin, size=Nsrcs
    ) ** -1
    spectral_indices = stats.norm.rvs(spectral_index_mean, spectral_index_std, Nsrcs)
    pix_indices = list(range(pix.size))
    intensities = np.zeros((pix.size, freqs.size), dtype=dtype)
    for intensity, index in zip(intensity_scales, spectral_indices):
        pixel = pix_indices.pop(stats.randint.rvs(0, len(pix_indices)))
        intensities[pixel,:] = intensity * (freqs / mfreq) ** index

    return intensities


def noiselike_eor(pix, freqs, amp=1e-5, seed=None, dtype=np.float32):
    """
    Mock up noiselike EoR intensities.

    Parameters
    ----------
    pix: np.ndarray
        Array with the same shape as the grid of sky pixels.
    freqs: np.ndarray of float
        Frequencies at which to evaluate the EoR signal.
    amp: float, optional
        Upper bound of the EoR intensity.
    seed: int, optional
        The random seed.
    dtype: type, optional
        Data type to use for simulating the foreground intensities

    Returns
    -------
    intensities: np.ndarray of float
        Shape ``(pix.size, freqs.size)`` array of intensities.
    """
    if seed is not None:
        np.random.seed(seed)
    return stats.uniform.rvs(0, amp, size=(np.size(pix), np.size(freqs)))
