#$ -S /bin/bash
#$ -N selfcal
#$ -j y
#$ -o grid_output/
#$ -V
#$ -cwd
#$ -l h_vmem=4G
#$ -l h_stack=256
ARGS=`pull_args.py $*`
CAL=psa455_v004_gc
NAME=PSA_gridtest
SCCAT=PAPER_pipeline_sky_model_v1
FLUXCAT=PAPER_pipeline_flux_cal_v1
FIMSTART=120
FIMSTOP=180
FIMDELT=6 #10 channels total
unset DISPLAY
echo casapy --logfile grid_output/casalog_selfcal_${JOB_ID}_${SGE_TASK_ID}.log \
--nologger -c ~/casascripts/bash_selfcal_test.py -C ${CAL} -s all \
--cat=${SCCAT} --fcat=${FLUXCAT} --prefix=${NAME} \
--fconfig=${FIMSTART}_${FIMSTOP}_${FIMDELT} ${ARGS}
casapy --logfile grid_output/casalog_selfcal_${JOB_ID}_${SGE_TASK_ID}.log \
--nologger -c  ~/casascripts/bash_selfcal_test.py -C ${CAL} -s all \
--cat=${SCCAT} --fcat=${FLUXCAT} --prefix=${NAME} \
--fconfig=${FIMSTART}_${FIMSTOP}_${FIMDELT}   ${ARGS}


