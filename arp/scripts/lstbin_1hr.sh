#! /bin/bash
LSTS=`python -c "import numpy as n; hr = n.pi/12; print ' '.join(['%f_%f' % (d,d+hr) for d in n.arange(0,24*hr,hr)])"`
echo $LSTS
for LST in $LSTS; do
    echo Working on $LST
    #lstbin.py -a 0_16 -p xx -C psa898_v002 --lst_res=42.95 --lst_rng=$LST --tfile=3600 $*
    lstbin.py -a cross -C psa898_v002 --lst_res=42.95 --lst_rng=$LST --tfile=600 $*
done;
