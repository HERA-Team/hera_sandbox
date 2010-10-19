#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output/
#$ -j y
#$ -N mdlvis
CAL=psa331_v005_gc
ARGS=`pull_args.py $*`

CATS=sumss_top1000
#SRCS=`srclist.py -s 10/0.15 --ra=9_21 --cat=${CATS} --divstr=, -C ${CAL}`

ssh djm growlnotify -m "Shredder_mesg  ${JOB_NAME}:${JOB_ID} "

echo mpirun -n $NSLOTS ~/scripts/mdl_vis_mpi_v002.py --nchan=1024 --sfreq=0.1 --sdf=0.00009 --startjd=2455335.2 --endjd=2455335.4 --pol=xx -C psa331_v005_gc -s cen,sgr,vir,hyd,her --inttime=5
mpirun -n $NSLOTS ~/scripts/mdl_vis_mpi_v002.py --nchan=1024 --sfreq=0.1 --sdf=0.00009 --startjd=2455335.2 --endjd=2455335.4 --pol=xx -C psa331_v005_gc -s cen,sgr,vir,hyd,her --inttime=5

ssh djm growlnotify -m "Shredder_mesg  ${JOB_NAME}:${JOB_ID} is finished"