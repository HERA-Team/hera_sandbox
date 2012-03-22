#$ -S /bin/bash
echo --daemons=\`qstat | qstat_to_hostport.py ${JOB_ID}\`
#---------------------------------------------------------------------------
fitmdl_daemon.py -p `python -c "print 53000+${SGE_TASK_ID}"`
