#! /bin/bash

FILE=$1

POL_ARR=5,6,7,8,11,12,18,19,24,25,27,28,29,30,32,33,34,35,36,39,48,51,52,55,60,75,76,77,78,79,91,92,93,94,95,107,108,109,110,111
HERA_HEX=9,10,20,22,31,43,53,64,65,72,80,81,88,89,96,97,104,105,112
PAPER_HEX=0,2,14,17,21,40,44,45,54,62,68,69,84,85,86,100,101,102,113
IMG_ARR=1,3,4,13,15,16,23,26,37,38,41,42,46,47,49,50,56,57,58,59,61,63,66,67,70,71,73,74,82,83,87,90,98,99,103,106,114,115,116,117,118,119,120,121,122,123,124,125,126,127

PATHY=$(dirname $FILE) #path to datafile
NAME=$(basename $FILE uv) #name of file without 'uv' at the end
BASE=${PATHY}/${NAME}
POL=${NAME:18:2}

echo 'working on ' ${FILE}

if [[ ! -e ${BASE}PP.uv ]]; then
    echo pull_antpols.py -p ${POL} -a "($POL_ARR)_($POL_ARR)" ${FILE}
    pull_antpols.py -p ${POL} -a "($POL_ARR)_($POL_ARR)" ${FILE}
    echo ${FILE}A -\> ${BASE}PP.uv
    mv ${FILE}A ${BASE}PP.uv
else
    echo ${BASE}PP.uv exists... skipping... 
fi 
   
if [[ ! -e ${BASE}HH.uv ]]; then 
    echo pull_antpols.py -p ${POL} -a "($HERA_HEX)_($HERA_HEX)" ${FILE}
    pull_antpols.py -p ${POL} -a "($HERA_HEX)_($HERA_HEX)" ${FILE}
    echo ${FILE}A -\> ${BASE}HH.uv
    mv ${FILE}A ${BASE}HH.uv
else
    echo ${BASE}HH.uv exists... skipping...
fi

if [[ ! -e ${BASE}PH.uv ]]; then
    echo pull_antpols.py -p ${POL} -a "($PAPER_HEX)_($PAPER_HEX)" ${FILE}
    pull_antpols.py -p ${POL} -a "($PAPER_HEX)_($PAPER_HEX)" ${FILE}
    echo ${FILE}A -\> ${BASE}PH.uv
    mv ${FILE}A ${BASE}PH.uv
else
    echo ${BASE}PH.uv exists... skipping...
fi

if [[ ! -e ${BASE}PI.uv ]]; then
    echo pull_antpols.py -p ${POL} -a "($IMG_ARR)_($IMG_ARR)" ${FILE}
    pull_antpols.py -p ${POL} -a "($IMG_ARR)_($IMG_ARR)" ${FILE}
    echo ${FILE}A -\> ${BASE}PI.uv
    mv ${FILE}A ${BASE}PI.uv
else
    echo ${BASE}PI.uv exists... skipping...
fi

    
    
