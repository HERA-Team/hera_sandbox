#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N ADD_FQ

f0=0.129354
df=0.00189
CAL=psa819_v007

for i in `seq 0 31`; do
    #fq=`python -c "print $f0+$i*$df"`
    #for pol in xx xy yx yy; do
    #    cl_img.py --tol=1e-6 -r radial --minuv=20 `ls -d snaps/zen*FB${i}${pol}_0000.dim.fits` 
    #    lin2stokes.py -C ${CAL} -f ${fq} \
    #        `ls -d snaps/zen*FB${i}${pol}_0000.bim.fits`
    #done
    #QUrot.py `ls -d snaps/zen*FB${i}[Q,U]_0000.bim.fits`
    for pol in I Qrot Urot V; do
        #mk_snap_map.py -f ${fq} -C ${CAL} -p xy \
        #    -i --nside=1024 -m synth_maps/RotSyn_FB${i}_${pol}0000.fits \
        #    `ls -d snaps/zen*FB${i}${pol}_0000.bim.fits`
        hpm_extract_facet.py -C ${CAL} -s "23:00_-30" --size=45_45 --res=4 -i \
            synth_maps/RotSyn_FB${i}_${pol}0000.fits
    done
done
