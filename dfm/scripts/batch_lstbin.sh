#! /bin/bash

LSTS='0_24'
FLEN='0.5'
lst_res=34
calfile=psa898_v003
files=$*
myCAPO=/home/damo/src/capo/dfm/scripts/
lst_binner=/home/damo/src/capo/zsa/scripts/lstbin_v02.py


for lbin in `${myCAPO}/gen_lsts_by_jid.py --lst_rng=${LSTS} --file_len=${FLEN}`; do
    echo python ${lst_binner} -C ${calfile} --lst_rng=${lbin} --lst_res=${lst_res} \
        --median --nsig=3 ${files}
    python ${lst_binner} -C ${calfile} --lst_rng=${lbin} --lst_res=${lst_res} \
        --median --nsig=3 ${files}
done
