### My Data ###
PREFIX='/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v3_xtalk/lstbin/'
EVEN_FILES=${PREFIX}'even/sep0,2/*I.uvGAL'
ODD_FILES=${PREFIX}'odd/sep0,2/*I.uvGAL'
RA='4_10' #'1_10'
NDAYS=20
EVEN_FILES=`lst_select.py -C psa6622_v003 --ra=${RA} ${EVEN_FILES[@]}`
ODD_FILES=`lst_select.py -C psa6622_v003 --ra=${RA} ${ODD_FILES[@]}`
SEP='0,2'
CALFILE='psa6622_v003'
CHAN='79_99'
NBOOT=20
POL='I'
DIRNAME='pspec_'${CHAN}'_'${SEP}'_'${RA}'_'${POL}

FRFEOR='--frfeor'

### Noise case instead of data ###
NOISE='' #'--noise_only'
NOISETYPE='' #'--same' #'--diff'
FRF='' #'--frf' 

### Other Options ###
MODE='' #'--oldsigloss'
LMODE='' #'--lmode=18'

#---------------------------------------------------------------------------------

mkdir ${DIRNAME}
echo Making Directory ${DIRNAME}
cd ${DIRNAME}

#Pspec cov
~/capo/pspec_pipeline/pspec_cov_v004.py --window=none -a cross -p I -c ${CHAN} -C ${CALFILE} -b ${NBOOT} ${NOISE} ${NOISETYPE} ${FRF} ${LMODE} ${EVEN_FILES} ${ODD_FILES}

#Bootstrap     
~/capo/pspec_pipeline/pspec_cov_boot_v002.py --identity pspec_boot*npz
mv pspec.npz pspec_I.npz
~/capo/pspec_pipeline/pspec_cov_boot_v002.py pspec_boot*npz
mv pspec.npz pspec_C.npz

#Sigloss
for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-1,3,10)))"` ; do
    mkdir inject_sep${SEP}_${inject}
    echo SIGNAL_LEVEL=${inject}
### OLD SIGLOSS ### 
    #~/capo/pspec_pipeline/sigloss_sim.py --window=none -a cross -p I -c ${CHAN} -C ${CALFILE} -b ${NBOOT} -i ${inject} ${NOISE} ${EVEN_FILES} ${ODD_FILES}
### NEW SIGLOSS ###
    ~/capo/pspec_pipeline/pspec_cov_v004_sigloss.py --window=none -a cross -p I -c ${CHAN} -C ${CALFILE} -b ${NBOOT} -i ${inject} ${NOISE} ${NOISETYPE} ${FRF} ${FRFEOR} ${MODE} ${LMODE} ${EVEN_FILES} ${ODD_FILES}
    echo ~/capo/pspec_pipeline/pspec_cov_v004_sigloss.py --window=none -a cross -p I -c ${CHAN} -C ${CALFILE} -b ${NBOOT} -i ${inject} ${NOISE} ${NOISETYPE} ${FRF} ${FRFEOR} ${MODE} ${LMODE} ${EVEN_FILES} ${ODD_FILES} > inject_sep${SEP}_${inject}/notes.txt
    mv *bootsigloss*.npz inject_sep${SEP}_${inject}/.
done

#Plotting & Noise Curve
~/capo/pspec_pipeline/plot_sigloss.py --output sigloss_factor
echo ~/capo/pspec_pipeline/plot_sigloss.py --output sigloss_factor
~/src/21cmSense/mk_array_file.py -C psa128 --bl_min=29.9 --bl_max=30.1 #XXX baselines hard-coded
echo ~/src/21cmSense/mk_array_file.py -C psa128 --bl_min=29.9 --bl_max=30.1 #XXX baselines hard-coded
~/src/21cmSense/calc_sense.py --ndays ${NDAYS} --n_per_day $[${RA##*_}-${RA%_*}]  --nchan 21 --bwidth 0.01 -m pess psa128.drift*npz
echo ~/src/21cmSense/calc_sense.py --ndays ${NDAYS} --n_per_day $[${RA##*_}-${RA%_*}]  --nchan 21 --bwidth 0.01 -m pess psa128.drift*npz
~/capo/pspec_pipeline/pspec_plot_simple.py --cov sigloss_factor.npz --SENSE *pess*.npz pspec_C.npz
echo ~/capo/pspec_pipeline/pspec_plot_simple.py --cov sigloss_factor.npz --SENSE *pess*.npz pspec_C.npz




