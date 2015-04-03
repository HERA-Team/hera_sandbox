#! /bin/bash
LSTS=`python -c "import numpy as n; hr = n.pi/12; print ' '.join(['%f_%f' % (d,d+hr) for d in n.arange(0,24*hr,hr)])"`
echo $LSTS
for LST in $LSTS; do
    echo Working on $LST
    lstbin.py -a cross -C psa898_v003 -s Sun --lst_res=42.95 --lst_rng=$LST --tfile=600 $*
done;
