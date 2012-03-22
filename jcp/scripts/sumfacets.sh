#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N sum_facets
#$ -o grid_output/

FACETS=`pull_facets.py --match="z.*xCd_(.*)_0000.dim.fits" $*`
FACET=`pull_args.py ${FACETS}`
echo $FACETS
echo $FACET
sum_fits.py z*uvcbrmMdxCd_${FACET}*.dim.fits -o facet_${FACET}.dim.fits
sum_fits.py z*uvcbrmMdxCd_${FACET}*.dbm.fits -o facet_${FACET}.dbm.fits