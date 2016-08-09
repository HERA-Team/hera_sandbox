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
