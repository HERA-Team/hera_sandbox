PREFIX='../../../lstbin_psa64_noise_1.5Jy_nofrf'
# PREFIX='../../lstbin_psa64_data_frf0'
# PREFIX='../../lstbin_psa64_data_frwidth_5_maxfr_5'
# PREFIX='../../lstbin_psa64_noise_3Jy_new_mdlvis_optimal'
EVEN_FILES=${PREFIX}'/even/sep0,1/lst*242.[3456]*.uvGAs'
ODD_FILES=${PREFIX}'/odd/sep0,1/lst*243.[3456]*.uvGAs'
# EVEN_FILES=${PREFIX}'/even/sep0,1/lst.24562*.[3456]*.uvGAL'
# ODD_FILES=${PREFIX}'/odd/sep0,1/lst.24562*.[3456]*.uvGAL'
# lst.24562*.[3456]
WD=$PWD #get the working directory
noise=''
boot=60
chans='95_115'
# rmbls='15_16,0_26,0_44,16_62,3_10,3_25'
# rmbls='7_12,7_52,8_11,11_36,14_40,14_54,15_16,16_62'
rmbls=''
t_eff=69
bl_length=30
# 95_115'
#export chans='30_50  51_71 78_98 95_115 103_123 127_147'
#chans='30_50 51_71 78_98 95_115 103_123 127_147'
#EVEN_FILES=`lst_select.py -C psa6622_v003 --ra=4_10 ${EVEN_FILES[@]}`
#ODD_FILES=`lst_select.py -C psa6622_v003 --ra=4_10 ${ODD_FILES[@]}`
SEP='0,1'
for chan in $chans; do
    # if [ $chan == 30_50 ]; then
    # continue
    # fi

    test -e $WD/${chan} || mkdir $WD/${chan}
    cd $WD/${chan}
    for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-2,-1,1)))"` ; do
        test -e inject_sep${SEP}_${inject} || mkdir inject_sep${SEP}_${inject}
        echo SIGNAL_LEVEL=${inject}

        #~/capo/pspec_pipeline/pspec_cov_v003_sigloss.py --window=none -a cross -p I -c 110_130 -C psa6622_v003 -b 20 -i ${inject} ${EVEN_FILES} ${ODD_FILES}
        #echo "~/capo/pspec_pipeline/pspec_cov_v003_sigloss.py --window=none -a cross -p I -c 110_130 -C psa6622_v003 -b 20 -i ${inject}" ${EVEN_FILES} ${ODD_FILES} > inject_sep${SEP}_${inject}/notes.txt

        #noise only
        /home/mkolopanis/src/capo/pspec_pipeline/sigloss_sim_no_gps.py --window=none -a cross -p I -c ${chan} -C psa6240_v003 -b ${boot} -i ${inject} ${noise} --rmbls=${rmbls} ${EVEN_FILES} ${ODD_FILES}
        echo "~/capo/pspec_pipeline/sigloss_sim_no_gps.py --window=none -a cross -p I -c ${chan} -C psa6240_v003 -b ${boot} -i ${inject} ${noise} --rmbls=${rmbls} "  ${EVEN_FILES} ${ODD_FILES} > inject_sep${SEP}_${inject}/notes.txt

        mv *bootsigloss*.npz inject_sep${SEP}_${inject}/.
    done
    cd $WD
done
for chan in $chans; do
    # continue
    cd $WD/${chan}
    python /home/mkolopanis/src/capo/pspec_pipeline/boots_to_pspec.py --t_eff=${t_eff} --bl_length=${bl_length}
    python /home/mkolopanis/src/capo/pspec_pipeline/sigloss_limits.py inject_sep*/pspec_pk_k3pk*.npz
    #cp sigloss.png ../sigloss_${chan}.png
    cd $WD #return to the sigloss dir to do the next channel
done
~/src/capo/mjk/scripts/plot_upper_lims_simple.py */pspec_limits_k3pk_p[CI]_85.npz --noisefiles='/home/mkolopanis/psa64/21cmsense_noise/dish_size_1/*drift_mod*.npz'   --outfile='noise_curve_85'
