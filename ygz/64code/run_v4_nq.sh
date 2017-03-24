#!bin/bash
python pspec_cov_v004.py --window=none --sep=sep0,1 --sepd=sep1,1 -b 20 -c 95_115 -p I -C psa6240_v003
python pspec_cov_v004.py --window=none --sep=sep0,1 --sepd=sep-1,1 -b 20 -c 95_115 -p I -C psa6240_v003
mv boot/* boot_nq/
python pspec_boot_multipath.py boot_nq/*
python plot_pk_k3pk.py out/pspec.npz
