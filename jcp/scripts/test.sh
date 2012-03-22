#$ -S /bin/bash
source .bashrc
#echo ${SGE_TASK_ID}
#echo ${JOB_ID}
#echo --daemons=\`qstat | qstat_to_hostport.py ${JOB_ID}\`
echo --daemons=`qstat | qstat_to_hostport.py ${JOB_ID}`
