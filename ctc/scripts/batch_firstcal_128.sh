#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=64G
#$ -l paper
#$ -N FIRSTCAL_128
ARGS=`pull_args.py $*`

echo ${ARGS}

for f in ${ARGS}; do
    if ((${f:33:7} < 2456679 )); then
        echo working on ${f}, which is in S1E1...
        #~/capo/omni/firstcal.py -C psa6622_v003 -p xx --outpath /data4/paper/2013EoR/Analysis/ProcessedData/epoch1/omni_v5_xtalk/ --ubls=64_10,64_49,64_3,64_41,64_25,64_19,64_48,64_29,64_24,64_28,64_55,64_34,64_27,64_51,64_57,64_9,64_22,64_20,64_43,64_53,64_31,10_65,10_72,10_80,10_88,10_96,10_104 --ex_ants=2,10,15,22,31,33,42,43,47,58,64,72,91,97,105,107,8,16,34,84,85,100 ${f} #S1E1xx
        #~/capo/omni/firstcal.py -C psa6622_v003 -p yy --outpath /data4/paper/2013EoR/Analysis/ProcessedData/epoch1/omni_v5_xtalk/ --ubls=64_10,64_49,64_3,64_41,64_25,64_19,64_48,64_29,64_24,64_28,64_55,64_34,64_27,64_51,64_57,64_9,64_22,64_20,64_43,64_53,64_31,10_65,10_72,10_80,10_88,10_96,10_104 --ex_ants=2,10,15,22,31,33,42,43,47,58,64,72,91,97,105,107,100,56,7,84,16 ${f} #S1E1yy
    fi
    if ((${f:33:7} > 2456678 )); then
        echo working on ${f}, which is in S1E2...
        #~/capo/omni/firstcal.py -C psa6622_v003 -p xx --outpath /data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v5_xtalk/ --ubls=64_10,64_49,64_3,64_41,64_25,64_19,64_48,64_29,64_24,64_28,64_55,64_34,64_27,64_51,64_57,64_9,64_22,64_20,64_43,64_53,64_31,10_65,10_72,10_80,10_88,10_96,10_104 --ex_ants=8,16,34,84,85,100 ${f} #S1E2xx
        #~/capo/omni/firstcal.py -C psa6622_v003 -p yy --outpath /data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v5_xtalk/ --ubls=64_10,64_49,64_3,64_41,64_25,64_19,64_48,64_29,64_24,64_28,64_9,64_22,64_20,64_43,64_53,64_31,10_65,10_72,10_80,10_88,10_96,10_104 --ex_ants=7,16,56,84,100 ${f} #S1E2yy # XXX had to eliminate some ubls to not get memory errors!
    fi
    if (( ${f:26:7} > 2456881 && ${f:26:7} < 2456929 )); then
        echo working on ${f}, which is in S2E3...
        #~/capo/omni/firstcal.py -C psa6622_v003 -p xx --outpath /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v5_xtalk/ --ubls=64_10,64_49,64_3,64_41,64_25,64_19,64_48,64_29,64_24,64_28,64_55,64_34,64_27,64_51,64_57,64_9,64_22,64_20,64_43,64_53,64_31,10_65,10_72,10_80,10_88,10_96,10_104 --ex_ants=8,13,16,24,26,34,37,38,85,107 ${f} #S2E3xx
        ~/capo/omni/firstcal.py -C psa6622_v003 -p yy --outpath /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v5_xtalk/ --ubls=64_10,64_49,64_3,64_41,64_25,64_19,64_48,64_29,64_24,64_28,64_55,64_34,64_27,64_51,64_57,64_9,64_22,64_20,64_43,64_53,64_31,10_65,10_72,10_80,10_88,10_96,10_104 --ex_ants=3,7,15,16,17,26,34,56,81,82,107 ${f} #S2E3yy
    fi
done
