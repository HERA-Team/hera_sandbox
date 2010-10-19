#! /bin/bash
DEST=/scratch/paper/pgb050/
for ((n=1; n <= 16; n++)); do echo "rsync -zav $* ${DEST}" | qsub -N rsync`python -c "print \"%02i\"%($n,)"` -j y -o grid_output/ -V -cwd -l h=node`python -c "print \"%02i\"%($n,)"`; done
