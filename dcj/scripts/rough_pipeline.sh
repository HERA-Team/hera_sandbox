#! /bin/bash
CAL=pgb966_v003
SUBSRCS=cyg,cas,crab,vir
DLYWID=5
DRTWID=5
CLN=1e-3
XRFI_CH=0_59,780_1023
NCHAN=256
AIPY=/usr/global/paper/bin/
# ---------------------------------------------------------
echo "using AIPY install in: "$AIPY
for F in $*; do
    echo Working on $F
    if ! ls ${F}dddddrxm &> /dev/null; then
        echo "Building ${F}ddrxa"
        ${AIPY}mdlvis.py -C $CAL -s $SUBSRCS -m sim $F
        ${AIPY}uv_addsub.py --sub $F ${F}s
        if ! ls ${F}drx.e_Sun &> /dev/null; then
            ${AIPY}xrfi.py -c $XRFI_CH -m val ${F}d
            ${AIPY}xtalk3.py ${F}dr
            ${AIPY}filter_src.py -C $CAL -e -p -s Sun -d $DLYWID -r $DRTWID --clean=$CLN ${F}drx
        fi
        if ! ls ${F}drx.e_cyg &> /dev/null; then
            ${AIPY}filter_src.py -C $CAL -e -p -s cyg -d $DLYWID -r $DRTWID --clean=$CLN ${F}drx
        fi        
        if ! ls ${F}drx.e_vir &> /dev/null; then
            ${AIPY}filter_src.py -C $CAL -e -p -s vir -d $DLYWID -r $DRTWID --clean=$CLN ${F}drx
        fi
        if ! ls ${F}drx.e_crab &> /dev/null; then
            ${AIPY}filter_src.py -C $CAL -e -p -s crab -d $DLYWID -r $DRTWID --clean=$CLN ${F}drx
        fi
        ${AIPY}uv_addsub.py --sub ${F}d ${F}drx.e_Sun
        ${AIPY}uv_addsub.py --sub ${F}dd ${F}drx.e_cyg
        ${AIPY}uv_addsub.py --sub ${F}ddd ${F}drx.e_vir
        ${AIPY}uv_addsub.py --sub ${F}dddd ${F}drx.e_crab
        ${AIPY}xrfi.py -c $XRFI_CH -m val ${F}ddddd
        ${AIPY}xtalk3.py ${F}dddddr
        #${AIPY}uv_addsub.py ${F}ddrx ${F}s
        combine_freqs.py -n $NCHAN ${F}dddddrx
        rm -rf ${F}d ${F}d{r,rx} ${F}dd ${F}d 
    fi
done
