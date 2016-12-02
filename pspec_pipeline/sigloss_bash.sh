### My Data ###
PREFIX='../..'
EVEN_FILES=${PREFIX}'/even/sep0,2/*I.uvGAL'
ODD_FILES=${PREFIX}'/odd/sep0,2/*I.uvGAL'
EVEN_FILES=`lst_select.py -C psa6622_v003 --ra=4_10 ${EVEN_FILES[@]}`
ODD_FILES=`lst_select.py -C psa6622_v003 --ra=4_10 ${ODD_FILES[@]}`
SEP='0,2'
CALFILE='psa6622_v003'
CHAN='110_130'
NBOOT=20

NOISE='' #'--noise_only'
NOISETYPE='' #'--same' #'--diff'
FRF='' #'--frf' 

FRFEOR='--frfeor'

MODE='' #'--oldsigloss'
LMODE='--lmode=17'

### Matt's Data
#EVEN_FILES='even/*uvGAL'
#ODD_FILES='odd/*uvGAL'
#SEP='0,1'
#CALFILE='psa6240_v003'
#CHAN='95_115'

#---------------------------------------------------------------------------------

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
    #~/capo/pspec_pipeline/pspec_cov_v003_sigloss.py --window=none -a cross -p I -c ${CHAN} -C ${CALFILE} -b ${NBOOT} -i ${inject} ${NOISE} ${EVEN_FILES} ${ODD_FILES}

    mv *bootsigloss*.npz inject_sep${SEP}_${inject}/.
done

#Plotting
~/capo/pspec_pipeline/plot_sigloss.py 



