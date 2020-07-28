"""Module for constructing visibilities/visibility fields."""

import numpy as np
from astropy import constants

from . import utils

def integrate_fringe_1d(pix, u_modes, dtype=np.complex64):
    """
    Integrate the fringe over each pixel, rather than sample it.
    
    Parameters
    ----------
    pix: array-like of float
        Coordinates of sky pixels partitioned in an equal-area way, with
        the zero point corresponding to zenith. Dimensionless.
    u_modes: array-like of float
        u-modes at which the sky is sampled.
    dtype: type, optional
        Data type to use when returning the integrated fringe.

    Returns
    -------
    fringe: np.ndarray of complex
        Fringe integrated over each pixel, with shape
        ``(pix.size, u_modes.size)``.
    """
    pix = np.array(pix).astype(dtype)
    u_modes = np.array(u_modes).astype(dtype)
    dx = np.mean(np.diff(pix))
    fringe_1 = np.exp(-2j * np.pi * np.outer(pix + dx / 2, u_modes))
    fringe_2 = np.exp(-2j * np.pi * np.outer(pix - dx / 2, u_modes))
    fringe = (fringe_1 - fringe_2) / (-2j * np.pi * np.outer(dx, u_modes))
    return fringe


def calculate_visibility_1d(
    antpos,
    pix,
    freqs,
    intensity,
    beam_power,
    dtype=np.complex64
):
    """
    Integrate the fringed, beam-weighted intensity over the sky at a single time.

    We use the following measurement equation to evaluate the visibilities:

    :math:`V = \\int A_\\nu(\\Omega)I_\\nu(\\Omega)\\exp(-i2\\pi\\u.n)d\\Omega`,

    where :math:`\\u = \\b\\nu/c`, :math:`A_\\nu(\\Omega)` is the peak-normalized
    primary beam response, :math:`I_\\nu(\\Omega)` is the specific intensity on
    the sky, and :math:`n` is the unit-vector pointing in the direction of solid
    angle :math:`\\Omega`.


    Parameters
    ----------
    antpos: dict
        Dictionary mapping antenna numbers to ENU positions, in meters.
    pix: array-like of float
        Coordinates of sky pixels partitioned in an equal-area way, with
        the zero point corresponding to zenith. Dimensionless.
    freqs: array-like of float
        Frequency at which the beam and sky intensity have been evaluated.
    intensity: array-like of float
        Sky intensity in units of Jy/beam, as a function of angular position and
        frequency.
    beam_power: array-like of float
        Peak-normalized beam power, as a function of angular position and
        frequency. (i.e. strictly non-negative values; square of E-field, peak
        normalized.)
    dtype: type, optional
        Data type to use when returning the visibilities.

    Returns
    -------
    visibility: np.ndarray of complex
        Visibilities calculated via the measurement equation, with shape
        ``(freqs.size, len(antpos), len(antpos))``.
    """
    antenna_numbers = list(antpos.keys())
    baselines = utils.get_baselines(antpos)
    intensity = np.atleast_2d(intensity).reshape(pix.size, freqs.size)
    beam_power = np.atleast_2d(beam_power).reshape(pix.size, freqs.size)
    dx = np.mean(np.diff(pix))
    visibilities = np.zeros((freqs.size, len(antpos), len(antpos)), dtype=dtype)
    for (ai, aj), (blx, bly, blz) in baselines.items():
        i = antenna_numbers.index(ai)
        j = antenna_numbers.index(aj)
        u_modes = blx * freqs * 1e9 / constants.c.value
        fringe = integrate_fringe_1d(pix=pix, u_modes=u_modes, dtype=dtype)
        vis_integrand = beam_power * intensity * fringe * dx
        V_ij = np.sum(vis_integrand, axis=0)
        visibilities[:,i,j] = V_ij
        visibilities[:,j,i] = V_ij.conj()
    return visibilities.astype(dtype)

