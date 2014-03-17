#$ -S /bin/bash
#$ -V
#$ -j y
#$ -cwd

ARGS=`pull_args.py /home/obs/Scratch/zen.*.uv`
ENDPATH=/home/obs/data/

for FILE in $ARGS
do
    echo ------------------------------------------------------------
    echo --- Working on ${FILE} ---
    echo ------------------------------------------------------------
    
    TRIPLET=`python /home/obs/Share/redux/get_uv_neighbor.py ${FILE}`
    echo "===> correct and first pass XRFI"
    if [ `python -c "print len('${TRIPLET}'.split())"` -lt 3 ]
    then
        echo No adjacent files to use. Skipping...
        continue
    fi
    for NEIGHBOR in $TRIPLET
    do
        if [ ! -f ${NEIGHBOR}cR ]
        then
            echo correct_and_XRFI.py -a 1 -t 80 -n 5 --df=6 -c 0_65,377_388,510,770,840,852,913,921_922,932_934,942_1023 ${NEIGHBOR}
            time correct_and_XRFI.py -a 1 -t 80 -n 5 --df=6 -c 0_65,377_388,510,770,840,852,913,921_922,932_934,942_1023 ${NEIGHBOR}
        fi
    done

    FILE=${FILE}cR
    TRIPLET=`python -c "t = '${TRIPLET}'.split(' '); print ' '.join([x+'cR' for x in t])"`

    echo "===> Improving RFI Flags"
    if [ !  -f ${FILE}R ]
    then
        #echo ddr_filter_coarse5.py -a 1 -p xx,xy,yx,yy --clean=1e-3 --maxbl=300 --output=ddr --invert $TRIPLET 
        #time ddr_filter_coarse5.py -a 1 -p xx,xy,yx,yy --clean=1e-3 --maxbl=300 --output=ddr --invert $TRIPLET 
        echo python /home/obs/Share/redux/ddr_filter_coarse.py -a 1 -p xx,xy,yx,yy --clean=1e-3 --maxbl=300 --output=ddr --invert $TRIPLET 
        time python /home/obs/Share/redux/ddr_filter_coarse.py -a 1 -p xx,xy,yx,yy --clean=1e-3 --maxbl=300 --output=ddr --invert $TRIPLET 
        echo xrfi_simple.py -a 1 --combine -t 40 -n 5 ${FILE}E --to_npz=${FILE}E.npz
        time xrfi_simple.py -a 1 --combine -t 40 -n 5 ${FILE}E --to_npz=${FILE}E.npz
        echo xrfi_simple.py -a all --combine -t 80 ${FILE} --from_npz=${FILE}E.npz
        time xrfi_simple.py -a all --combine -t 80 ${FILE} --from_npz=${FILE}E.npz
    fi
   
    TRIPLET=`python -c "t = '${TRIPLET}'.split(' '); print ' '.join([t[0], t[1]+'R', t[2]])"`
    echo "===> Beginning Compression"
    #echo ddr_filter_coarse5.py -a all -p xx,xy,yx,yy --clean=1e-3 --maxbl=300 --nsections=10 $TRIPLET
    #time ddr_filter_coarse5.py -a all -p xx,xy,yx,yy --clean=1e-3 --maxbl=300 --nsections=10 $TRIPLET
    echo python /home/obs/Share/redux/ddr_filter_coarse.py -a all -p xx,xy,yx,yy --clean=1e-3 --maxbl=300 --nsections=20 $TRIPLET
    time python /home/obs/Share/redux/ddr_filter_coarse.py -a all -p xx,xy,yx,yy --clean=1e-3 --maxbl=300 --nsections=20 $TRIPLET

    echo "===> Moving data to its final resting place"
    echo scp -rp -c arcfour256 ${FILE}R[DEF] qmaster:${ENDPATH}
    time scp -rp -c arcfour256 ${FILE}R[DEF] qmaster:${ENDPATH}
    echo scp -rp -c arcfour256 ${FILE}E.npz qmaster:${ENDPATH}
    time scp -rp -c arcfour256 ${FILE}E.npz qmaster:${ENDPATH}
done
