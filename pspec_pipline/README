Stuff from Delay Spec paper:
gen_psrc_spec_npz.py:
    - generate mock npz foregrounds used for plot_k3pk_vs_k_vs_fq
gen_sync_spec_npz.py:
    - ditto
plot_k3pk_vs_k_vs_fq.py:
    - Run the delay spectrum on mock npz data
    - Plot dspec vs freq waterfall
plot_k3pk_vs_kperp_vs_kpara.py:
    - Run the delay spectrum on mock npz data
    - plot dspec vs kpara vs kperp for a particular freq

Actual data pipeline stuff:
pspec_prep.py:
    - do the full-band fg clean, needs to be updated to k3pk approach
    - should probably save teh clean components to a parallel uv file
    - add uv variable 'bin' for each baseline/time, so it doesn't have to be computed
    - might convert to one half of the uv plane here & conjugate if necessary
    - should change to use the capo.dspec module
pspec_to_npz_quick_n_dirty.py
    - read uv files and bin together measurements, save to npz files
    - should modify to use uv['bin'] from pspec_prep.py
    - accumulates to separate "slots" in each bin for cross-multiplying w/o bias
plot_pspec_npz.py:
    - reads in npz files from pspec_to_npz_quick_n_dirty
    - runs remaining sub-band delay spectrum stuff
    - squares measurements within bins
    - accumulates all bin measurements into a final k3pk power spectrum plot
pspec_pipeline.sh:
    - recommended way to reduce data from uvcbt files
pspec_npz_bls_in_bin.py:
    - allows you to query specific bin numbers in npz files generated above.


Old stuff that's probably obselete:
uv2pk.py
