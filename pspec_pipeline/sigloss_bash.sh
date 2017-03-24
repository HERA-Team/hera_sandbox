<<<<<<< HEAD
PREFIX='../../lstbin_psa64_data_frf0'
EVEN_FILES=${PREFIX}'/even/sep0,1/lst*242.[3456]*.uvGAL'
ODD_FILES=${PREFIX}'/odd/sep0,1/lst*243.[3456]*.uvGAL'
WD=$PWD #get the working directory
noise=''
boot=60
chans='30_50 95_115'
# 95_115'
#export chans='30_50  51_71 78_98 95_115 103_123 127_147'
#chans='30_50 51_71 78_98 95_115 103_123 127_147'
#EVEN_FILES=`lst_select.py -C psa6622_v003 --ra=4_10 ${EVEN_FILES[@]}`
#ODD_FILES=`lst_select.py -C psa6622_v003 --ra=4_10 ${ODD_FILES[@]}`
SEP='0,1'
for chan in $chans; do
    continue
    test -e $WD/${chan} || mkdir $WD/${chan}
    cd $WD/${chan}
    for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-1,2,50)))"` ; do
        mkdir inject_sep${SEP}_${inject}
        echo SIGNAL_LEVEL=${inject}

        #~/capo/pspec_pipeline/pspec_cov_v003_sigloss.py --window=none -a cross -p I -c 110_130 -C psa6622_v003 -b 20 -i ${inject} ${EVEN_FILES} ${ODD_FILES}
        #echo "~/capo/pspec_pipeline/pspec_cov_v003_sigloss.py --window=none -a cross -p I -c 110_130 -C psa6622_v003 -b 20 -i ${inject}" ${EVEN_FILES} ${ODD_FILES} > inject_sep${SEP}_${inject}/notes.txt

        #noise only
        /home/mkolopanis/src/capo/pspec_pipeline/sigloss_sim.py --window=none -a cross -p I -c ${chan} -C psa6240_v003 -b ${boot} -i ${inject} ${noise}  ${EVEN_FILES} ${ODD_FILES}
        echo "~/capo/pspec_pipeline/sigloss_sim.py --window=none -a cross -p I -c ${chan} -C psa6240_v003 -b ${boot} -i ${inject} ${noise} " ${EVEN_FILES} ${ODD_FILES} > inject_sep${SEP}_${inject}/notes.txt

        mv *bootsigloss*.npz inject_sep${SEP}_${inject}/.
    done
    cd $WD
done
for chan in $chans; do
    cd $WD/${chan}
    /home/mkolopanis/src/capo/pspec_pipeline/plot_sigloss_boots.py
    cp sigloss.png ../sigloss_${chan}.png
    cd $WD #return to the sigloss dir to do the next channel
done
=======
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




>>>>>>> 1bf4764e9898309bbe32ec5d3c34f3f87da9aa38
