# PREFIX='../../lstbin_psa64_noise_1.5Jy_optimal'
PREFIX='lstbin_psa64_data_frf0'
# PREFIX='../../lstbin_psa64_data_frwidth_5_maxfr_5'
# PREFIX='../../lstbin_psa64_noise_3Jy_new_mdlvis_optimal'
# EVEN_FILES=${PREFIX}'/even/sep0,1/lst*242.[3456]*.uvGAL'
# ODD_FILES=${PREFIX}'/odd/sep0,1/lst*243.[3456]*.uvGAL'
PREFIX='../../../'$PREFIX
EVEN_FILES=${PREFIX}'/even/sep0,1/lst.24562*.[3456]*.uvGAL'
ODD_FILES=${PREFIX}'/odd/sep0,1/lst.24562*.[3456]*.uvGAL'
# lst.24562*.[3456]
WD=$PWD #get the working directory
RUN_NAME='sigloss_sim_data_jacknife_test'
noise=''
boot=60
chans='30_50 95_115'
rmbls_list='15_16,0_26,0_44,16_62,3_10,3_25'
bls_master='1_4 1_58 2_33 2_43 4_17 5_18 5_32 6_33 6_52 7_12 7_52 8_11 8_45 9_58 11_36 12_38 13_17 13_56 14_40 14_54 15_21 18_35 20_63 21_53 22_61 23_30 24_48 24_55 25_48 27_55 27_57 28_29 28_34 30_32 31_45 34_51 35_61 36_60 39_46 39_60 41_47 41_49 42_63 44_62 56_59'

t_eff=69
bl_length=30
# 95_115'
#export chans='30_50  51_71 78_98 95_115 103_123 127_147'
#chans='30_50 51_71 78_98 95_115 103_123 127_147'
#EVEN_FILES=`lst_select.py -C psa6622_v003 --ra=4_10 ${EVEN_FILES[@]}`
#ODD_FILES=`lst_select.py -C psa6622_v003 --ra=4_10 ${ODD_FILES[@]}`
SEP='0,1'


run_dir=$WD/$RUN_NAME
test -e $run_dir || mkdir $run_dir
cd  $run_dir

echo -e 'This is sigloss_bash_v3.sh'

for bl in ${bls_master}; do
    bl_dir=$run_dir/rmbl_${bl}
    test -e $bl_dir || mkdir $bl_dir
    cd $bl_dir
    rmbls=$rmbls_list,$bl

    echo -e 'Removing Baselines ' $rmbls

    for chan in $chans; do
        # continue
        chan_dir=$bl_dir/$chan
        test -e $chan_dir || mkdir $chan_dir
        cd ${chan_dir}

        echo -e '\tChannels ' $chan

        for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-3,4,60)))"` ; do
            inj_dir=$chan_dir/inject_sep${SEP}_${inject}

            test -e $inj_dir || mkdir $inj_dir
            echo SIGNAL_LEVEL=${inject}
            LOGFILE=$inj_dir/sigloss.log



            /home/mkolopanis/src/capo/pspec_pipeline/sigloss_sim.py --window=none -a cross -p I -c ${chan} -C psa6240_v003 -b ${boot} -i ${inject} ${noise} --rmbls=${rmbls} ${EVEN_FILES} ${ODD_FILES} | tee -a $LOGFILE
            echo "~/capo/pspec_pipeline/sigloss_sim.py --window=none -a cross -p I -c ${chan} -C psa6240_v003 -b ${boot} -i ${inject} ${noise} --rmbls=${rmbls} "  ${EVEN_FILES} ${ODD_FILES} > $inj_dir/notes.txt

            mv *bootsigloss*.npz $inj_dir
        done
        cd  $bl_dir
    done
    for chan in $chans; do
        cd  $bl_dir/${chan}
        python /home/mkolopanis/src/capo/pspec_pipeline/boots_to_pspec.py --t_eff=${t_eff} --bl_length=${bl_length}
        python /home/mkolopanis/src/capo/pspec_pipeline/sigloss_limits.py inject_sep*/pspec_pk_k3pk*.npz
        #cp sigloss.png ../sigloss_${chan}.png
        cd  $bl_dir #return to the sigloss dir to do the next channel
    done

    ~/src/capo/mjk/scripts/plot_noise_curves.py */pspec_limits_k3pk_p[CI]_85.npz --noisefiles='../../../21cmsense_noise/new_noise/psa6240_v003.drift_mod*.npz'   --outfile='noise_curve_rmbl_'$bl'_85'

done
