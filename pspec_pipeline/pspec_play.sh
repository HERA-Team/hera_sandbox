PREFIX='..'
EVEN_FILES=${PREFIX}'/even/sep0,2/*I.uvGAL'
ODD_FILES=${PREFIX}'/odd/sep0,2/*I.uvGAL'
#EVEN_FILES=${PREFIX}'/even/sep0,2/*I.uvGA'
#ODD_FILES=${PREFIX}'/odd/sep0,2/*I.uvGA'

EVEN_FILES=`lst_select.py -C psa6622_v003 --ra=4_10 ${EVEN_FILES[@]}`
ODD_FILES=`lst_select.py -C psa6622_v003 --ra=4_10 ${ODD_FILES[@]}`
SEP='0,2'
CALFILE='psa6622_v003'
CHAN='110_130'
NBOOT=20

TRICK='C_original_longtime'
#TRICK='C_Cnofrf'

#CHANGE=''
CHANGE='--changeC'

#CHANGETYPE=''
#CHANGETYPE='--reg=100' #requires --changeC
#CHANGETYPE='--otherbls="0,2"' #requires --changeC
#CHANGETYPE='--CnoFRF #requires --changeC
#CHANGETYPE='--Cfg' #requires --changeC
CHANGETYPE='--Clongtime' #requires --changeC, can also be used in combination with --Cfg or --CnoFRF

FRFEOR='--frfeor'

#--------------------------------------------------------------------------------

#Pspec & Signal Loss
for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-1,3,10)))"` ; do
    mkdir inject_sep${SEP}_${inject}
    echo SIGNAL_LEVEL=${inject}
    ~/capo/pspec_pipeline/pspec_cov_v004_play.py --window=none -a cross -p I -c ${CHAN} -C ${CALFILE} -b ${NBOOT} ${CHANGE} ${CHANGETYPE} ${FRFEOR} -i ${inject} ${EVEN_FILES} ${ODD_FILES}
    mv *bootsigloss*.npz inject_sep${SEP}_${inject}/.
done

#Bootstrap     
~/capo/pspec_pipeline/pspec_cov_boot_v002.py --identity pspec_boot*npz
mv pspec.npz pspec_I.npz
~/capo/pspec_pipeline/pspec_cov_boot_v002.py pspec_boot*npz
mv pspec.npz pspec_C.npz

#Plotting
~/capo/pspec_pipeline/plot_sigloss.py --output 'factor'

#Move files into directory
rm -r ${TRICK}
mkdir ${TRICK}
mv inject* pspec* sigloss* factor* ${TRICK}
