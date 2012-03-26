#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N sum_facets
#$ -o grid_output/
#$ -l h_vmem=0.5G

FACETS=`pull_facets.py $* --match="z.*_(.\d+)_0000.dim.fits"`
FACETS=`pull_args.py ${FACETS}`
echo $FACETS
for FACET in $FACETS
do
echo $FACET
print_pointing.py z*_${FACET}_0000.dim.fits
sum_fits.py z*_${FACET}_0000.dim.fits -o facet_${FACET}.dim.fits
sum_fits.py z*_${FACET}_0000.dbm.fits -o facet_${FACET}.dbm.fits
done