#$ -S /usr/bin/bash
#$ -N phsgrid
#$ -j y

SRC=cen
CAL=psa331_v004_gc
mpirun -n $NSLOTS ~/scripts/phsgrid_mpi_2.py -pxx -s ${SRC} -C ${CAL} $* --exp