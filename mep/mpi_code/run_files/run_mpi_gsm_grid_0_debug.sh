#PBS -q debug
#PBS -l mppwidth=24
#PBS -l walltime=00:30:00
#PBS -N mc_0_db
#PBS -e out_files/mc_0_db_.$PBS_JOBID.err
#PBS -o out_files/mc_0_db_.$PBS_JOBID.out
#PBS -V

echo "before module loads"

module load python
module load mpi4py

echo "after module loads"

cd $PBS_O_WORKDIR 
# num0,beam_sig,del_bl,num_bl (for one side of grid)
aprun -n 24 python-mpi /global/homes/m/mpresley/scripts/monte_carlo/mpi_monte_carlo_gen_y.py 5 0.087 4 10

echo "aprun was submitted"
