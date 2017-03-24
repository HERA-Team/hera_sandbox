#!bin/bash
python pspec_cov_v004.py --data_dir=/data2/ali_et_al_2015_apj_data/ --window=none --sep=sep0,1 --sepd=sep0,1 -b 20 -c 95_115 -p I -C psa6240_v003
python pspec_cov_v004.py --data_dir=/data2/ali_et_al_2015_apj_data/ --window=none --sep=sep-1,1 --sepd=sep-1,1 -b 20 -c 95_115 -p I -C psa6240_v003
python pspec_cov_v004.py --data_dir=/data2/ali_et_al_2015_apj_data/ --window=none --sep=sep1,1 --sepd=sep1,1 -b 20 -c 95_115 -p I -C psa6240_v003
python pspec_boot_v2.py boot/sep0,1_sep0,1/* boot/sep1,1_sep1,1/* boot/sep-1,1_sep-1,1/*
python plot_pk_k3pk.py out/pspec.npz
