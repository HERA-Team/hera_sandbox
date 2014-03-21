#$ -S /bin/bash
#$ -V 
#$ -cwd
#$ -j y
#$ -N pspec
#$ -o grid_output/

dir=$(pull_args.py $1)
echo "Working on ${dir}"
cd $dir
_chn=$(echo ${dir} | awk '{split($1,a,"/"); print a[2]}')
_pol=$(echo ${dir} | awk '{split($1,a,"/"); print a[3]}')
_sep=$(echo ${dir} | awk '{split($1,a,"/"); print a[4]}')
${SCRIPTSDIR}/pspec_redmult_cov_gps.py -b ${NBOOT} -a ${_sep} -c ${_chn} -p ${_pol} ${FILES} \ 
    && ${SCRIPTSDIR}/pspec_pk_k3pk_boot.py pspec_boot*npz 
cd -
