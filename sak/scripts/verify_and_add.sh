#!/bin/bash

source /usr/global/paper/CanopyVirtualEnvs/PAPER_Distiller/bin/activate

for FILE in $*; do
    python check_corruption.py $FILE  && echo $FILE
done >> MYGOODFILES.txt
add_observations.py `cat MYGOODFILES.txt`
