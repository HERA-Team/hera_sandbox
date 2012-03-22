#! /bin/bash
CAL=pgb966_v003
SUBSRCS=cyg,cas,crab,vir
DLYWID=5
DRTWID=5
CLN=1e-3
XRFI_CH=0_59,780_1023

# ---------------------------------------------------------

for F in $*; do
    echo Working on $F
    if ! ls ${F}ddrxa &> /dev/null; then
        echo "Building ${F}ddrxa"
        mdlvis.py -C $CAL -s $SUBSRCS --sim $F
        uv_addsub.py --sub $F ${F}s
        if ! ls ${F}drx.e_Sun &> /dev/null; then
            xrfi.py -c $XRFI_CH -m val ${F}d
            xtalk3.py ${F}dr
            filter_src.py -C $CAL -e -p -s Sun -d $DLYWID -r $DRTWID --clean=$CLN ${F}drx
        fi
        uv_addsub.py --sub ${F}d ${F}drx.e_Sun
        xrfi.py -c $XRFI_CH -m val ${F}dd
        xtalk3.py ${F}ddr
        uv_addsub.py ${F}ddrx ${F}s
        rm -rf ${F}d ${F}d{r,rx} ${F}dd ${F}dd{r,rx}
    fi
done
