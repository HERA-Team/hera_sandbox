#! /bin/bash
#launches a single task server a node
#$ -V
#$ -cwd
#$ -l paper
#$ -l h_vmem=15G
#$ -j y
#$ -o grid_output
mkdir /scratch/${SGE_TASK_ID}
still_taskserver.py /scratch/${SGE_TASK_ID}
