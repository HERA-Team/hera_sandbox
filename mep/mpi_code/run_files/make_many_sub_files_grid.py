import sys
import os

# kk=0
# for beam_sig in (0.087,0.175,0.349,0.689,1.047):
#     for del_bl in (8,10,20):#4,6)
#         fcontent = """#PBS -q regular
# #PBS -l mppwidth=48
# #PBS -l walltime=05:00:00
# #PBS -N Q_grid_{0}
# #PBS -e out_files/Q_grid_{1}.$PBS_JOBID.err
# #PBS -o out_files/Q_grid_{2}.$PBS_JOBID.out
# #PBS -V

# module load python
# module load mpi4py

# cd $PBS_O_WORKDIR 
# # beam_sig,del_bl,num_bl (for one side of grid)
# aprun -n 48 python-mpi /global/homes/m/mpresley/scripts/mpi_Q_matrix_grid.py {3} {4} 10
# """.format(kk,kk,kk,beam_sig,del_bl)
#         with open('./run_mpi_Q_matrix_grid_{0}.sh'.format(kk), 'w') as file:
#             file.writelines(fcontent)
#         os.system('qsub ./run_mpi_Q_matrix_grid_{0}.sh'.format(kk))
#         kk += 1

kk=0
# for beam_sig in (0.087,0.175,0.349,0.689,1.047):
#     for del_bl in (8,10,20):#4,6)
for beam_sig, del_bl in ((0.69,20),(1.05,20),(0.69,10),(1.05,10),(0.69,8),(1.05,8),(0.09,6),(0.09,4)):
    fcontent = """#PBS -q regular
#PBS -l mppwidth=72
#PBS -l walltime=05:00:00
#PBS -N gsm_grid_{0}
#PBS -e out_files/gsm_grid_{1}.$PBS_JOBID.err
#PBS -o out_files/gsm_grid_{2}.$PBS_JOBID.out
#PBS -V

module load python
module load mpi4py

cd $PBS_O_WORKDIR 
# beam_sig,del_bl,num_bl (for one side of grid)
aprun -n 72 python-mpi /global/homes/m/mpresley/scripts/mpi_gsm_grid.py {3} {4} 10
""".format(kk,kk,kk,beam_sig,del_bl)
    with open('./run_mpi_gsm_grid_{0}.sh'.format(kk), 'w') as file:
        file.writelines(fcontent)
    os.system('qsub ./run_mpi_gsm_grid_{0}.sh'.format(kk))
    kk += 1
