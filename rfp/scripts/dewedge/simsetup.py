"""Tools for setting up a simulation."""

import yaml

import astropy_healpix as aph
import numpy as np

from . import antpos
from . import beams
from . import sky

def parse_config(config_file):
    """Parse a configuration file and do minimal processing."""
    with open(config_file, "r") as cfg:
        config = yaml.load(cfg.read(), Loader=yaml.FullLoader)

    # Set up the antenna array.
    array_params = config["setup"]["array_layout"]
    array_type = array_params.pop("array_type", "golomb")
    array_layout = getattr(antpos, f"{array_type}_array")(**array_params)
    
    # Set up the frequency and sky coordinate axes.
    freq_params = config["setup"]["freq"]
    freqs = np.linspace(
        freq_params["fmin"], freq_params["fmax"], freq_params["Nfreqs"]
    ).astype(freq_params["dtype"])
    sky_params = config["setup"]["sky"]
    pix = sky_coords_from_params(**sky_params)
    
    # Construct the beams.
    beam_params = config["beams"]
    beam_dict = {}
    if type(pix) is np.ndarray:
        ndim = 1
        beam = beams.Beam1d(xs=pix, freqs=freqs)
    else:
        ndim = 2
        beam = beams.Beam2d(pix=pix, freqs=freqs)
    for beam_key, beam_kwargs in beam_params.items():
        beam_type = beam_kwargs.pop("type")
        beam_dict[beam_key] = getattr(beam, beam_type)(**beam_kwargs)

    # Construct the intensities.
    sky_params = config["sky"]
    intensities = {}
    for source, kwargs in sky_params.items():
        key = kwargs.pop("label", source)
        simulator = getattr(sky, source)
        intensities[key] = simulator(pix=pix, freqs=freqs, **kwargs)
        
    # Collect objects and return.
    metadata = {
        "antpos": array_layout,
        "freqs": freqs,
        "pix": pix,
        "ndim": ndim,
        "vis_dtype": config["setup"].get("vis_dtype", np.complex64),
    }
    parsed_config = {
        "metadata": metadata,
        "filing": config["filing"],
        "beams": beam_dict,
        "intensities": intensities,
    }
    return parsed_config


def sky_coords_from_params(ndim, **kwargs):
    """Create an array of sky coordinates."""
    if ndim == 1:
        sky_coords = np.linspace(
            kwargs.get("xmin", -1), kwargs.get("xmax", 1), kwargs.get("Npix", 500),
        ).astype(kwargs.get("dtype", np.float32))
    elif ndim == 2:
        pass  # placeholder
    else:
        raise NotImplementedError("Only 1-d skies currently supported.")

    return sky_coords
