#! /bin/bash

for d in `ls -d psa*`
do
    echo $d `ls -1d ${d}/*uv | wc -l`
done
