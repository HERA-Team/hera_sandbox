#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=8G

ARGS=`pull_args.py /home/obs/Scratch/zen.*.uvcR`
for FILE in $ARGS; do
    TRIPLET=`python /home/obs/Share/redux/get_uv_neighbor.py ${FILE}`
    if [ `echo ${TRIPLET} | wc -w` != 3 ]; then
        echo "No adjacent files to use. Skipping..."
        continue
    fi
    echo "===> Improving RFI Flags"
    if [ ! -f ${FILE}R ]; then
        echo python /home/obs/Share/redux/ddr_filter_coarse.py -a 1 -p xx,xy,yx,yy --clean=1e-3 --maxbl=300 --output=ddr --invert $TRIPLET 
        time python /home/obs/Share/redux/ddr_filter_coarse.py -a 1 -p xx,xy,yx,yy --clean=1e-3 --maxbl=300 --output=ddr --invert $TRIPLET
        echo xrfi_simple.py -a 1 --combine -t 80 -n 5 ${FILE}E --to_npz=${FILE}E.npz
        time xrfi_simple.py -a 1 --combine -t 80 -n 5 ${FILE}E --to_npz=${FILE}E.npz
        echo xrfi_simple.py -a all --combine -t 80 ${FILE} --from_npz=${FILE}E.npz
        time xrfi_simple.py -a all --combine -t 80 ${FILE} --from_npz=${FILE}E.npz
    fi
done
