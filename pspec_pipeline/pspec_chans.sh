#! /bin/bash
PREFIX=Feb7_2
mkdir $PREFIX
#chans=`python -c "print ' '.join(['%d_%d'%(i,i+39) for i in range(10,150,1)])"`
chans='37_70 80_119 100_139 110_149 119_150'
RA=1:01_9:00
NBOOT=20
#DATAPATH=fringe_hor_v006
SCRIPTSDIR=~/scripts/
DATAPATH=`pwd`
PIDS=""
for chan in $chans;
do
if [ ! -e ${PREFIX}/${chan}/pspec_${PREFIX}_${chan}.png ]
then
mkdir ${PREFIX}/${chan}

#BEGIN EW BLS
echo Starting EW
mkdir ${PREFIX}/${chan}/pspec_EW
cd ${PREFIX}/${chan}/pspec_EW
${SCRIPTSDIR}/pspec_redmult_cov_gps.py -b $NBOOT -a 0_16  -c $chan \
 `~/scripts/lst_select.py -C psa898_v003 --ra=${RA} ${DATAPATH}/lst*uv` && ~/scripts/pspec_pk_k3pk_boot.py \
 pspec_boot*npz &
PIDS="${PIDS} "$!
cd -

mkdir ${PREFIX}/${chan}/pspec_NE
cd ${PREFIX}/${chan}/pspec_NE
#BEGIN NE BLS
echo Starting NE
${SCRIPTSDIR}/pspec_redmult_cov_gps.py -b $NBOOT -a 1_16 -c $chan \
 `~/scripts/lst_select.py -C psa898_v003 --ra=${RA} ${DATAPATH}/lst*uv` && ~/scripts/pspec_pk_k3pk_boot.py \
 pspec_boot*npz &
PIDS="${PIDS} "$!
cd - 


mkdir ${PREFIX}/${chan}/pspec_SE
cd ${PREFIX}/${chan}/pspec_SE
echo Starting SE
${SCRIPTSDIR}/pspec_redmult_cov_gps.py -b $NBOOT -a 0_17 -c $chan \
 `~/scripts/lst_select.py -C psa898_v003 --ra=${RA} ${DATAPATH}/lst*uv` && ~/scripts/pspec_pk_k3pk_boot.py \
 pspec_boot*npz &
PIDS="${PIDS} "$!
cd -

fi

done
echo waiting on power spectra ${PIDS}
wait $PIDS

echo averaging power spectra for channels
for chan in $chans;
do
echo $chan
#PLOT
/home/jacobsda/scripts/pspec_plot_pk_k3pk.py ${PREFIX}/${chan}/pspec_??/pspec.npz
mv pspec_pk_k3pk.npz ${PREFIX}/${chan}/
mv pspec.png pspec_${PREFIX}_${chan}.png
cp  pspec_${PREFIX}_${chan}.png ${PREFIX}/${chan}/
done

