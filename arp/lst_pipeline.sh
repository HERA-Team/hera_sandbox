#! /bin/bash
CAL=pgb015_v004
SUBSRCS=cyg,cas,crab,vir
DLYWID=5
DRTWID=5
CLN=1e-3
XRFI_CH=0_240,780_1023

# ---------------------------------------------------------

for F in $*; do
    echo Working on $F
    if ! ls ${F}ddrxda &> /dev/null; then
        echo "Building ${F}ddrxa"
        mdlvis.py -C $CAL -s $SUBSRCS $F
        uv_addsub.py --sub $F ${F}s
        mdlvis.py -C $CAL -s Sun ${F}d
        uv_addsub.py --sub ${F}d ${F}ds
        xrfi.py -c $XRFI_CH -m val ${F}dd
        xtalk3.py ${F}ddr
        filter_src.py -C $CAL -e -p -s Sun -d $DLYWID -r $DRTWID --clean=$CLN ${F}ddrx
        uv_addsub.py --sub ${F}ddrx ${F}ddrx.e_Sun
        uv_addsub.py ${F}ddrx.e_Sun ${F}ds
        uv_addsub.py ${F}ddrxd ${F}s
        rm -rf ${F}{s,d} ${F}d{s,d} ${F}dd{r,rx,rxd} ${F}ddrx.e_Sun
    fi
done
