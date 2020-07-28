"""Script for experimenting with foreground subtraction."""

import argparse
import datetime
import logging
import os
import sys
sys.path.append("/home/bobby/HERA/hera_sandbox/rfp/scripts/")

import numpy as np
from astropy import constants
from scipy.interpolate import interp1d

from uvtools.dspec import dayenu_filter

from dewedge import simsetup
from dewedge import utils
from dewedge import vis

parser = argparse.ArgumentParser(
    description="Mock up some visibilities and attempt to subtract off the foregrounds."
)
parser.add_argument("config", type=str, help="Path to configuration file.")
parser.add_argument("-v", "--verbose", action="store_true", help="Print updates.")
parser.add_argument("--clobber", default=False, action="store_true", help="Clobber files.")
args = parser.parse_args()

if __name__ == "__main__":
    logger = logging.getLogger("fg_sub_script")

    # Setup.
    config = simsetup.parse_config(args.config)
    antpos = config["metadata"]["antpos"]
    freqs = config["metadata"]["freqs"]
    pix = config["metadata"]["pix"]
    ndim = config["metadata"]["ndim"]
    vis_dtype = config["metadata"]["vis_dtype"]
    beams = config["beams"]
    intensities = config["intensities"]
    filing_info = config["filing"]
    log_filename = filing_info.get("log_filename", "fg_sub.log")
    outdir = filing_info.get("outdir", ".")
    metadata_filename = filing_info.get("metadata_filename", "metadata")
    metadata_filename = os.path.join(outdir, metadata_filename)
    if args.verbose:
        logging.basicConfig(
            filename=os.path.join(outdir, log_filename),
            filemode="a",
            format="%(message)s",
            level=logging.DEBUG,
        )

    # Write some useful information to disk now, so we can delete it to save space.
    # Beam information:
    beam_filename = filing_info.get("beam_filename", "beams")
    beam_filename = os.path.join(outdir, beam_filename)
    if utils.save_file(beam_filename, args.clobber):
        np.savez(beam_filename, **beams)

    # Intensities:
    intensity_filename = filing_info.get("intensities_filename", "intensities")
    intensity_filename = os.path.join(outdir, intensity_filename)
    if utils.save_file(intensity_filename, args.clobber):
        np.savez(intensity_filename, **intensities)
    
    # Construct visibilities.
    logger.info("Constructing visibilities...")
    visibilities = {}
    simulator = getattr(vis, f"calculate_visibility_{ndim}d")
    for beam_key, beam_power in beams.items():
        for sky_key, intensity in intensities.items():
            full_key = beam_key + sky_key
            logger.info(f"Making visibilities of type {full_key}...")
            logger.info(f"Starting at {datetime.datetime.now()}")
            logger.info(".............................................")
            visibilities[full_key] = simulator(
                antpos=antpos,
                pix=pix,
                freqs=freqs,
                intensity=intensity,
                beam_power=beam_power,
                dtype=vis_dtype,
            )
            logger.info(f"Finished at {datetime.datetime.now()}")
            logger.info(".............................................")

    logger.info("Generating summed visibilities...")
    logger.info(f"Starting at {datetime.datetime.now()}")
    logger.info(".............................................")
    summed_visibilities = {}
    for beam_key in beams:
        for sky_key in intensities:
            if sky_key == "E":
                continue
            eor_key = beam_key + "E"
            vis_key = beam_key + sky_key
            sum_key = vis_key + "E"
            summed_visibilities[sum_key] = visibilities[vis_key] + visibilities[eor_key]
    logger.info(f"Finished at {datetime.datetime.now()}")
    del beams, intensities

    # Write the true visibilities to disk; delete them to save space.
    logger.info("Writing true visibilities to disk and removing them from memory...")
    vis_filename = filing_info.get("vis_filename", "true_visibilities")
    vis_filename = os.path.join(outdir, vis_filename)
    if utils.save_file(vis_filename, args.clobber):
        np.savez(vis_filename, **visibilities)
    logger.info(f"True visibilities written to disk at {vis_filename}")
    del visibilities, config

    # Do foreground subtraction on each of the summed visibilities.
    logger.info("Now performing foreground subtraction...")
    fg_subtracted_visibilities = {}
    antnums = list(antpos.keys())
    baselines = utils.get_baselines(antpos)
    baselines = {
        antpair: bl * 1e9 / constants.c.value for antpair, bl in baselines.items()
    }  # Convert to nanoseconds.
    flags = np.zeros_like(list(summed_visibilities.values())[0]).astype(np.bool)
    coverage = np.zeros_like(flags).astype(np.int16)
    minimum_coverage = 10  # Make this a parameter we can control.
    # TODO: add support for specifying filter parameters in a configuration file
    filter_center = 0
    filter_half_width = 40
    filter_factor = 1e-9
    filter_kwargs = {
        "filter_centers": [filter_center],
        "filter_half_widths": [filter_half_width],
        "filter_factors": [filter_factor],
        "filter_dimensions": [0],
    }

    # TODO: see if there's a faster way to do this
    for vis_key, visibilities in summed_visibilities.items():
        logger.info(f"Beginning foreground subtraction for visibility type {vis_key}")
        logger.info(f"Starting at {datetime.datetime.now()}")
        logger.info(".............................................")
        fg_sub_vis = np.zeros_like(visibilities)
        for (ant1, ant2), bl_vec in baselines.items():
            logger.info(f"Beginning foreground subtraction for antpair {(ant1,ant2)}")
            logger.info(f"Starting at time {datetime.datetime.now()}")
            logger.info(".............................................")
            i = antnums.index(ant1)
            j = antnums.index(ant2)
            for ind, freq in enumerate(freqs):
                uvw = freq * bl_vec
                active_freqs = []
                interp_vis = []
                for (ai, aj), (blx, bly, blz) in baselines.items():
                    # TODO: extend this to acting on uv-modes rather than just u-modes
                    afreq = uvw[0] / blx
                    if afreq < freqs.min() or afreq > freqs.max():
                        continue  # Don't have access to the information we need.
                    coverage[ind,i,j] += 1
                    active_freqs.append(afreq)
                    # Make the interpolator.
                    i_ = antnums.index(ai)
                    j_ = antnums.index(aj)
                    re_spline = interp1d(
                        freqs, visibilities[:,i_,j_].real, kind="cubic", assume_sorted=True
                    )
                    im_spline = interp1d(
                        freqs, visibilities[:,i_,j_].imag, kind="cubic", assume_sorted=True
                    )
                    interp_vis.append(re_spline(afreq) + 1j * im_spline(afreq))
                
                # Flag the channel if poorly constrained.
                well_covered = coverage[ind,i,j] >= minimum_coverage
                out_of_bounds = freq < min(active_freqs) or freq > max(active_freqs)
                if out_of_bounds or not well_covered:
                    flags[ind,i,j] = True
                    continue

                # Sort the arrays in case the baselines aren't ordered by length.
                sort_key = np.argsort(active_freqs)
                active_freqs = np.asarray(active_freqs)[sort_key]
                interp_vis = np.asarray(interp_vis)[sort_key]

                # Apply a low-pass filter to keep only the foregrounds at this u.
                foreground_vis = interp_vis - dayenu_filter(
                    x=active_freqs,
                    data=interp_vis,
                    wgts=np.ones_like(active_freqs),
                    **filter_kwargs
                )[0]

                # Subtract the foregrounds, then record the result for this baseline/freq.
                interp_vis -= foreground_vis
                re_spline = interp1d(
                    active_freqs, interp_vis.real, kind="cubic", assume_sorted=True
                )
                im_spline = interp1d(
                    active_freqs, interp_vis.imag, kind="cubic", assume_sorted=True
                )
                fg_sub_vis[ind,i,j] = re_spline(freq) + 1j * im_spline(freq)

        # Finally, record the result.
        fg_subtracted_visibilities[vis_key] = fg_sub_vis

    # Now write everything to disk.
    logger.info(f"Writing results to disk at {datetime.datetime.now()}")
    logger.info(".............................................")
    
    # Foreground-subtracted data:
    fg_sub_filename = "_".join([os.path.splitext(vis_filename)[0], "foreground_subtracted"])
    fg_sub_filename = fg_sub_filename + os.path.splitext(vis_filename)[1]
    if utils.save_file(fg_sub_filename, args.clobber):
        np.savez(fg_sub_filename, **fg_subtracted_visibilities)

    # Simulation/Analysis metadata:
    metadata = {}
    metadata["antenna_numbers"] = np.array(antnums)
    metadata["antenna_positions"] = np.array(list(antpos.values()))
    metadata["freqs"] = freqs
    metadata["pix"] = pix
    metadata["ndim"] = ndim
    metadata.update(
        {"flags": flags, "coverage": coverage, "minimum_coverage": minimum_coverage}
    )
    metadata.update(filter_kwargs)
    if utils.save_file(metadata_filename, args.clobber):
        np.savez(metadata_filename, **metadata)

