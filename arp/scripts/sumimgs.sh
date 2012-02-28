#! /usr/bin/env bash
#CH_LIST="240_360 360_480 480_600 600_720"
CH_LIST="240 360 480 600 719"
#TYPE_LIST="xRF RFF"
TYPE_LIST="RFF"
for CH in $CH_LIST; do
    for TYPE in $TYPE_LIST; do
        ../scripts/sumimg.py 1330_4000_01*${TYPE}_c${CH}.rim.fits
        mv sumimg.fits 1330_4000_sum${TYPE}_c${CH}.rim.fits
        ../scripts/sumimg.py 1330_4000_01*${TYPE}_c${CH}.dbm.fits
        mv sumimg.fits 1330_4000_sum${TYPE}_c${CH}.dbm.fits
    done
done
