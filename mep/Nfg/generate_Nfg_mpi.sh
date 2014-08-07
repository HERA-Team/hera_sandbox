# !/bin/bash
#PBS -S /bin/bash
#PBS -N generate_Nfg_mpi
#PBS -j eo
#PBS -l mppwidth=24,walltime=01:30:00
#PBS -q regular
#PBS -A m1871

module swap PrgEnv-pgi PrgEnv-gnu
module load gcc
module load openmpi
module load python/2.7.3
module load numpy
module load matplotlib
module list

#nside=8
nside=32
nlmax=200
ClFname="Cl_silly.dat"
#templateFname="template_small.dat"
templateFname="template.dat"
KfgFname="Kfg.npy"
GmatrixFname="Gmatrix.npy"
calFile="basic_amp_aa"
numAnts=8
freq=100 # in MHz
beamSig=0.3 # in radians
NfgFname="Nfg.npy"

cd $GSCRATCH/GlobalSignalInterferometer/Nfg

#gfortran generate_Kfg.f90 -o generate_Kfg.x -O3
#mpif90 generate_Kfg_mpi.f90 -o generate_Kfg_mpi.x -O3
#ftn -O3 -ffast-math -funroll-loops -o generate_Kfg_mpi.x generate_Kfg_mpi.f90


rm -f args.dat
touch args.dat
echo $nside >> args.dat
echo $nlmax >> args.dat
echo $ClFname >> args.dat
echo $templateFname >> args.dat
echo Kfg.raw >> args.dat

echo "Forming the foreground covariance matrix in image space..."
#./generate_Kfg.x args.dat
#mpirun -np 7 generate_Kfg_mpi.x args.dat
time aprun -n 24 -N 24 ./generate_Kfg_mpi.x args.dat
time python raw2npy.py Kfg.raw $KfgFname
echo "Forming the instrumental response matrices..."
time python generate_G.py $calFile $numAnts $nside $freq $beamSig $GmatrixFname
echo "Multiplying matrices together..."
time python sandwich.py $nside $KfgFname $GmatrixFname $NfgFname
echo "Done!"
