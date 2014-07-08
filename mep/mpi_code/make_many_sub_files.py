import sys
import os

kk=0
for beam_sig in (0.087,0.175,0.349,0.689,1.047):
    for del_bl in (4,6,8):#,10,20):
        fcontent = """#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=03:30:00
#PBS -N grid_{0}
#PBS -e grid_{1}.$PBS_JOBID.err
#PBS -o grid_{2}.$PBS_JOBID.out
#PBS -V

module load python
module load mpi4py

cd $PBS_O_WORKDIR 
# beam_sig,del_bl,num_bl (for one side of grid)
aprun -n 24 python-mpi /global/homes/m/mpresley/scripts/many_grid_runs/mpi_Q_matrix_grid.py {3} {4} 10
""".format(kk,kk,kk,beam_sig,del_bl)
        with open('./run_file/run_mpi_Q_matrix_{0}.sh'.format(kk), 'w') as file:
            file.writelines(fcontent)
        os.system('qsub ./run_file/run_mpi_Q_matrix_{0}.sh'.format(kk))
        kk += 1
