echo 'cleaning out scratch'
for n in 02 03 04 05 06 07 08 09 10 ; do ssh node${n} rm -rf /scratch/paper/*; done
echo 'du -csh /scratch/paper'
for n in 02 03 04 05 06 07 08 09 10; do ssh node${n} du -csh /scratch/paper; done

