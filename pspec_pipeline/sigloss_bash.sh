### My Data ###
PREFIX='..'
EVEN_FILES=${PREFIX}'/even/sep0,2/frf_new/*I.uvGAL'
ODD_FILES=${PREFIX}'/odd/sep0,2/frf_new/*I.uvGAL'
EVEN_FILES=`lst_select.py -C psa6622_v003 --ra=4_10 ${EVEN_FILES[@]}`
ODD_FILES=`lst_select.py -C psa6622_v003 --ra=4_10 ${ODD_FILES[@]}`
SEP='0,2'
CALFILE='psa6622_v003'
chan='110_130'
noise='--noise_only'

### Matt's Data
#EVEN_FILES='even/lst*242.[3456]*'
#ODD_FILES='odd/lst*243.[3456]*'
#SEP='0,1'
#CALFILE='psa6240_v003'
#chan = '30_50'
#noise = ''


for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-3,-1,5)))"` ; do
    mkdir inject_sep${SEP}_${inject}
    echo SIGNAL_LEVEL=${inject}
    
    ~/capo/pspec_pipeline/sigloss_sim.py --window=none -a cross -p I -c ${chan} -C ${CALFILE} -b 20 -i ${inject} ${noise} ${EVEN_FILES} ${ODD_FILES}
    echo "~/capo/pspec_pipeline/sigloss_sim.py --window=none -a cross -p I -c ${chan} -C ${CALFILE} -b 20 -i ${inject} ${noise}" ${EVEN_FILES} ${ODD_FILES} > inject_sep${SEP}_${inject}/notes.txt

    mv *bootsigloss*.npz inject_sep${SEP}_${inject}/.
done

~/capo/pspec_pipeline/plot_sigloss.py
