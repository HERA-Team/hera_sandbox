import sys
import os

# kk=0
# #for beam_sig in (0.087,0.175,0.349,0.689,1.047):
# for del_bl in (10,20,4,6,8):
#     fcontent = """#PBS -q regular
# #PBS -l mppwidth=24
# #PBS -l walltime=02:30:00
# #PBS -N Q_hy_{0}
# #PBS -e out_files/Q_hy_{1}.$PBS_JOBID.err
# #PBS -o out_files/Q_hy_{2}.$PBS_JOBID.out
# #PBS -V

# module load python
# module load mpi4py

# cd $PBS_O_WORKDIR 
# # del_bl,num_bl (for one side of grid)
# aprun -n 24 python-mpi /global/homes/m/mpresley/scripts/mpi_Q_matrix_hybrid_grid.py {3} 10
# """.format(kk,kk,kk,del_bl)
#     with open('./run_mpi_Q_matrix_{0}.sh'.format(kk), 'w') as file:
#         file.writelines(fcontent)
#     os.system('qsub ./run_mpi_Q_matrix_{0}.sh'.format(kk))
#     kk += 1

kk=0
for del_bl in (10,20,4,6,8):
    fcontent = """#PBS -q regular
#PBS -l mppwidth=72
#PBS -l walltime=05:00:00
#PBS -N gsm_hy_{0}
#PBS -e out_files/gsm_hy_{1}.$PBS_JOBID.err
#PBS -o out_files/gsm_hy_{2}.$PBS_JOBID.out
#PBS -V

module load python
module load mpi4py

cd $PBS_O_WORKDIR 
# del_bl,num_bl (for one side of grid)
aprun -n 72 python-mpi /global/homes/m/mpresley/scripts/mpi_gsm_hybrid.py {3} 10
""".format(kk,kk,kk,del_bl)
    with open('./run_mpi_gsm_{0}.sh'.format(kk), 'w') as file:
        file.writelines(fcontent)
    os.system('qsub ./run_mpi_gsm_{0}.sh'.format(kk))
    kk += 1
