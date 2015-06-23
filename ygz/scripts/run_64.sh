#!/bin/bash

scDIR=
dataDIR=/Users/yunfanzhang/local/simuDATA/pspec_2456249.24936.uv
cuefile=P0_15.cue
N=10
lsf1=$(ls $DIR | grep pspec_2)
lsf2=$(ls $DIR | grep pspec_0_38)


cd $scDIR
while read -r line
do
    name=$line
    echo "Name read from file - $name"
done < "$cuefile"

python ./64.py $filename $filename --ant 0_26_0_26 --chan 100 --pol xx --lst 2456249.278_2456249.278
