#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=12G
#$ -l paper
#$ -N POSTPROCESS_128
ARGS=`pull_args.py $*`

CALFILE='psa6622_v003'
echo ${ARGS}

for f in ${ARGS}; do
    ~/capo/omni/omni_xrfi.py --boxside=16 ${f:65:20}.npz 
    ~/capo/omni/omni_xrfi_apply.py --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v4_xtalk/%s.chisqflag.npz ${f}
    xrfi_simple.py -n 1000 -c 101,102,148,149,150,151,152,153,154,169 ${f}R
    ~/capo/ctc/scripts/pspec_prep.py -C ${CALFILE} -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${f}RR
    xrfi_simple.py -n 4 ${f}RRB 
done
