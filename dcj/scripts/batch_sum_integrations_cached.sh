#$ -S /bin/bash
cd /scratch/paper/psa112/
rm -r *uvcbI
sum_integrations.py -n 10 *uvcb
