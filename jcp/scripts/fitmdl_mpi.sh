#!/bin/sh
# OpenMPI example
# Export all environment variables
#$ -V
# Your job name
#$ -N test
# Use current working directory
#$ -cwd
# Join stdout and stder
#$ -j y
# PARALLEL ENVIRONMENT:
# -pe ompi 32
# Enable resource reservation
#$ -R y
# The max hard walltime for this job is 9 hours (after this it will be killed)
#$ -l h_rt=09:00:00
# The max soft walltime for this job is 10 hours (after this SIGUSR2 will be sent)
#$ -l s_rt=10:00:00
# The following is for reporting only. It is not really needed
# to run the job. It will show up in your output file.
echo "Got $NSLOTS processors."
# The mpirun command.
mpirun -np $NSLOTS /data1/paper/pober/scripts/fitmdl_mpi_v002.py  $*
#mpirun -np $NSLOTS /home/jacobsda/scripts//beamform_mpi_test.py -s J0950+142,J1002+285,J1115+404,J1143+221,J1145+313,J1423+194 --cat=helm_flagged --gsm=jys,index,dra,ddec -C pgb015_v005_gc -pxx --maxiter=500 --quiet -b -d 10 --clean=5e-4 --name=helm_flagged_20jy_10%ion --dt=4n --time helm_flagged_20jy_10%ion_2.uv
