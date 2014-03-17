#$ -S /bin/bash
#$ -V
#$ -l h_vmem=8G
#$ -cwd

ARGS=`pull_args.py /home/obs/Scratch/zen.*.uvcR`
ENDPATH=/home/obs/data/

for FILE in $ARGS; do
    TRIPLET=`python /home/obs/Share/redux/get_uv_neighbor.py ${FILE}`
    RTRIPLET=""
    for t in $TRIPLET; do
        [[ -e ${t}R ]] && RTRIPLET="${RTRIPLET} ${t}R" || RTRIPLET="${RTRIPLET} ${t}"
    done
    echo python /home/obs/Share/redux/ddr_filter_coarse.py -a all -p xx,xy,yx,yy --clean=1e-3 --maxbl=300 --nsections=20 ${RTRIPLET}
    time python /home/obs/Share/redux/ddr_filter_coarse.py -a all -p xx,xy,yx,yy --clean=1e-3 --maxbl=300 --nsections=20 ${RTRIPLET}

    echo "===> Moving data to its final resting place"
    echo scp -rp -c arcfour256 ${FILE}R[DEF] qmaster:${ENDPATH}
    time scp -rp -c arcfour256 ${FILE}R[DEF] qmaster:${ENDPATH}
    echo scp -rp -c arcfour256 ${FILE}E.npz qmaster:${ENDPATH}
    time scp -rp -c arcfour256 ${FILE}E.npz qmaster:${ENDPATH}
done
